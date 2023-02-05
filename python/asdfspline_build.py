from subprocess import run

from cffi import FFI

run(
    ['cargo', 'build', '--all', '--release'],
    cwd='..',
    check=True,
)

run(['cbindgen', '../ffi', '--config', 'cbindgen.toml', '-o', 'asdfspline.h'],
    check=True)

ffibuilder = FFI()
ffibuilder.cdef(open('asdfspline.h').read())
ffibuilder.set_source(
    '_asdfspline',
    r"""
    #include <stdbool.h>
    #include "asdfspline.h"
    """,
    include_dirs=['.'],
    libraries=['asdfspline_ffi'],
    library_dirs=['../target/release'],
)

if __name__ == '__main__':
    ffibuilder.compile(verbose=True)
