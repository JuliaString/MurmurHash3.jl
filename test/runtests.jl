# This file includes code that was formerly a part of Julia.
# License is MIT: LICENSE.md

@static VERSION < v"0.7.0-DEV" ? (using Base.Test) : (using Test)

@test mmhash128(SubString("--hello--",3,7)) == mmhash128("hello")
