[![Code Check](https://github.com/stefan-kurtz/gttl/actions/workflows/code-check.yml/badge.svg)](https://github.com/stefan-kurtz/gttl/actions/workflows/code-check.yml)
[![Linux CI](https://github.com/stefan-kurtz/gttl/actions/workflows/make-tests-linux.yml/badge.svg)](https://github.com/stefan-kurtz/gttl/actions/workflows/make-tests-linux.yml)
[![MacOS CI](https://github.com/stefan-kurtz/gttl/actions/workflows/make-tests-macos.yml/badge.svg)](https://github.com/stefan-kurtz/gttl/actions/workflows/make-tests-macos.yml)
[![Windows CI](https://github.com/stefan-kurtz/gttl/actions/workflows/make-tests-windows.yml/badge.svg)](https://github.com/stefan-kurtz/gttl/actions/workflows/make-tests-windows.yml)
# gttl
Sources of the genometools template library

Note that Windows users will need to install ZLib as an additional dependency, eg. by running
```
vcpkg install zlib
```
With a `vcpkg` instance in `C:\vcpkg`.
Windows users will also need to set the environment variable `CXX="clang++"`.
The build system will automatically find the headers and library files. If `vcpkg` is installed in a non-standard location, you can set the environment variable `VCPKG_ROOT` to your vcpkg installation path.
