from setuptools import setup

setup(
    name='asdfspline',
    version='0.0.0',
    author='Matthias Geier',
    author_email='Matthias.Geier@gmail.com',
    description='Python bindings for ASDF splines',
    long_description=open('README.rst').read(),
    package_dir={'': 'src'},
    py_modules=['asdfspline'],
    cffi_modules=['asdfspline_build.py:ffibuilder'],
    setup_requires=["CFFI>=1.0.0"],
    install_requires=["CFFI>=1.0.0"],
    zip_safe=False,
)
