# This file includes code that was formerly a part of Julia.
# License is MIT: LICENSE.md

using MurmurHash3

using Test

p1 = SubString("--hello--",3,7)
p2 = "hello"

_memhash(siz, ptr) = ccall(Base.memhash, UInt, (Ptr{UInt8}, Csize_t, UInt32), ptr, siz, 0%UInt32)
mh(str::String) = _memhash(sizeof(str), pointer(str))
mh(str::AbstractString) = mh(string(str))

mmhash(str::String) = mmhash128_a(sizeof(str), pointer(str), 0%UInt32)
mmhashc(str::AbstractString) = mmhash128_c(str, 0%UInt32)

mh32(str) = mmhash32(sizeof(str), pointer(str), 0%UInt32)

@testset "MurmurHash3" begin
    @test mmhashc(p1) == mmhash(p2)
    @static if sizeof(Int) == 8
        @test last(mmhashc(p1)) == mh(p1)
        @test last(mmhashc(p2)) == mh(p2)
        @test last(mmhash(p2))  == mh(p1)
    else
        @test mh32(p2) == mh(p2)
    end
end


