"""
mmap-based FITS output write — workaround attempt for the macOS U-state drizzle hang.

The hang is a kernel lost-wakeup in the buffered-write path: numpy's
``ndarray.tofile`` issues a ``write()`` whose ``copyin`` (under ``cluster_write``)
page-faults and parks on an rw-lock that is never woken (``lck_rw_sleep``). The
process wedges unkillably (~17-26 s CPU in) and only a reboot clears it. Chunking
the ``write()`` at the astropy layer did NOT help (2026-07-16) — the race is per
``write()`` syscall, in ``cluster_write``'s ``copyin``.

This routes the big data write around ``cluster_write``/``copyin`` entirely: it
``mmap``s the output file region and copies the array bytes in with a numpy store
(a ``memcpy``), so the write goes through the ``vm_fault``-on-mapped-file path
instead. The mmap *read* path already works fine on this machine, so the mmap
*store* path has a real chance of dodging the lost-wakeup.

``install()`` monkeypatches ``astropy.io.fits.file._array_to_file`` so every real
on-disk-file data write goes through the mmap path. No-op off macOS. Import and
call ``install()`` before AstroDrizzle runs; because the drizzle scripts re-exec
themselves for the no-CR subprocess pass, installing at import covers both.

Gotcha handled: astropy opens output files write-only, but ``mmap`` with write
access needs an ``O_RDWR`` fd, so we open a second ``O_RDWR`` fd on the same path
for the mapping and keep astropy's buffered stream consistent via ``flush()`` /
``seek()`` around the mmap write.
"""

import os
import sys
import mmap

import numpy as np
import astropy.io.fits.file as _fitsfile
from astropy.io.fits.util import isfile

# Copy the array into the mapping in bounded pieces so the memcpy touches the
# file-backed destination pages incrementally rather than in one huge fault burst.
_COPY_BYTES = int(os.environ.get('DRIZZLE_MMAP_COPY_MB', '8')) * 1024 * 1024
_orig_array_to_file = _fitsfile._array_to_file
_installed = False


def _mmap_array_to_file(arr, outfile):
    """mmap+memcpy replacement for astropy's ``_array_to_file`` for real files."""
    try:
        seekable = outfile.seekable()
    except AttributeError:
        seekable = False

    # Only real on-disk files hit the wedging cluster_write path and can be
    # mmap'd; anything else (BytesIO, pipes, non-seekable) uses astropy's original.
    name = getattr(outfile, 'name', None)
    if not (isfile(outfile) and seekable and isinstance(name, str)):
        return _orig_array_to_file(arr, outfile)

    # Contiguous raw-byte view of the array in its current (already big-endian,
    # on-disk) memory order — identical bytes to what arr.tofile(f) would write.
    src = np.ascontiguousarray(arr.view(np.ndarray)).reshape(-1).view(np.uint8)
    nbytes = int(src.size)
    if nbytes == 0:
        return

    # Flush astropy's buffered header bytes so the fd offset and file length are
    # truthful before we map, then record where the data must land.
    outfile.flush()
    pos = outfile.tell()

    # mmap offsets must be a multiple of the allocation granularity (page size);
    # map from the page-aligned start at/below pos and copy at the sub-page delta.
    pagesize = mmap.ALLOCATIONGRANULARITY
    map_start = (pos // pagesize) * pagesize
    delta = pos - map_start
    maplen = delta + nbytes

    # Separate O_RDWR fd so PROT_WRITE mmap works even though astropy opened the
    # file write-only. Same inode/page-cache as astropy's fd.
    fd = os.open(name, os.O_RDWR)
    try:
        if os.fstat(fd).st_size < pos + nbytes:
            os.ftruncate(fd, pos + nbytes)
        mm = mmap.mmap(fd, maplen, flags=mmap.MAP_SHARED,
                       prot=mmap.PROT_READ | mmap.PROT_WRITE, offset=map_start)
        try:
            dest = np.ndarray((nbytes,), dtype=np.uint8, buffer=mm, offset=delta)
            step = max(pagesize, _COPY_BYTES)
            for i in range(0, nbytes, step):
                dest[i:i + step] = src[i:i + step]
            mm.flush()
            if os.environ.get('DRIZZLE_MMAP_DEBUG'):
                print(f'[mmap_fits_write] mmap-wrote {nbytes} bytes to '
                      f'{os.path.basename(name)} at offset {pos}', file=sys.stderr, flush=True)
        finally:
            mm.close()
    finally:
        os.close(fd)

    # Reposition astropy's stream past the data so its trailing padding write
    # (a sub-2880-byte tail, harmless) continues at the right offset.
    outfile.seek(pos + nbytes)


def install():
    """Activate the mmap-write patch. macOS only; safe to call repeatedly."""
    global _installed
    if sys.platform != 'darwin':
        return False
    if not _installed:
        _fitsfile._array_to_file = _mmap_array_to_file
        _installed = True
    return True
