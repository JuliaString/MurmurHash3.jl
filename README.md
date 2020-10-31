# MurmurHash3: a pure Julia implementation of the MurmurHash3 functions

[pkg-url]: https://github.com/JuliaString/MurmurHash3.jl.git

[julia-url]:    https://github.com/JuliaLang/Julia
[julia-release]:https://img.shields.io/github/release/JuliaLang/julia.svg

[release]:      https://img.shields.io/github/release/JuliaString/MurmurHash3.jl.svg
[release-date]: https://img.shields.io/github/release-date/JuliaString/MurmurHash3.jl.svg

[license-img]:  http://img.shields.io/badge/license-MIT-brightgreen.svg?style=flat
[license-url]:  LICENSE.md

[gitter-img]:   https://badges.gitter.im/Join%20Chat.svg
[gitter-url]:   https://gitter.im/JuliaString/Lobby?utm_source=badge&utm_medium=badge&utm_campaign=pr-badge

[travis-url]:   https://travis-ci.org/JuliaString/MurmurHash3.jl
[travis-s-img]: https://travis-ci.org/JuliaString/MurmurHash3.jl.svg
[travis-m-img]: https://travis-ci.org/JuliaString/MurmurHash3.jl.svg?branch=master

[codecov-url]:  https://codecov.io/gh/JuliaString/MurmurHash3.jl
[codecov-img]:  https://codecov.io/gh/JuliaString/MurmurHash3.jl/branch/master/graph/badge.svg

[contrib]:    https://img.shields.io/badge/contributions-welcome-brightgreen.svg?style=flat

[![][release]][pkg-url] [![][release-date]][pkg-url] [![][license-img]][license-url] [![contributions welcome][contrib]](https://github.com/JuliaString/MurmurHash3.jl/issues)

| **Julia Version** | **Unit Tests** | **Coverage** |
|:------------------:|:------------------:|:---------------------:|
| [![][julia-release]][julia-url] | [![][travis-s-img]][travis-url] | [![][codecov-img]][codecov-url]
| Julia Latest | [![][travis-m-img]][travis-url] | [![][codecov-img]][codecov-url]

Note: this can be used to replace the C MurmurHash library used by base Julia,
to implement a `hash` function that gives compatible results.

This provides the following functions (written in pure Julia):
`mmhash128` which hashes a UTF-8 (or ASCII, which is compatible) string
`mmhash128_c` which can be used to hash an generic abstract string (by converting to UTF-8 on the fly, without having to allocate the string, working a few characters at a time)
Note that the hashes on 64-bit systems are *not* the same as on 32-bit systems
(this is true for base Julia `hash` as well)

`mmhash32` which creates a 32 bit hash from a UTF-8/ASCII string
(this works for `String`, as well as some of the `Str` types, such as `ASCIIStr`, `UTF8Str`, `Binary`, `Text1Str` can all be hashed directly)
(note, there is no mmhash32_c yet, so other strings have to be converted to `String` type before hashing)

Julia uses the following code to create a hash using the C MurmurHash library:
```
const memhash = UInt === UInt64 ? :memhash_seed : :memhash32_seed
const memhash_seed = UInt === UInt64 ? 0x71e729fd56419c81 : 0x56419c81

function hash(s::String, h::UInt)
    h += memhash_seed
    ccall(memhash, UInt, (Ptr{UInt8}, Csize_t, UInt32), s, sizeof(s), h % UInt32) + h
end
```
similar code, such as that in `JuliaString/StrBase.jl/src/hash.jl`, implements the
`hash` function for all of the `Str` types, but is specialized for performance based on the type of the string and whether or not it is aligned in memory.

