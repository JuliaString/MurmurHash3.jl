# This file includes code that was formerly a part of Julia.
# License is MIT: LICENSE.md

using MurmurHash3
using MurmurHash3: mmhash128_4, mmhash128_8_c, mmhash128_8_a, mmhash128_8_u

using Test

p1 = SubString("--hello--",3,7)
p2 = "hello"

const MaybeSub = Union{String,SubString{String}}
const mhseed32 = 0x56419c81
const mhseed64 = 0x71e729fd56419c81

# 32-bit MurmurHash (calling C version in Base)
mh32c(siz::Int, ptr::Ptr{UInt8}, seed=0%UInt32) =
    ccall(:memhash32_seed, UInt32, (Ptr{UInt8}, Csize_t, UInt32), ptr, siz, seed)
mh32c(str::MaybeSub, seed=0%UInt32) = mh32c(sizeof(str), pointer(str), seed)

# Calling Julia version in MurmurHash.jl
mh32j(str::MaybeSub, seed=0%UInt32) = mmhash32(sizeof(str), pointer(str), seed)
mh32j(len::Integer, val::Unsigned, seed=0%UInt32) = mmhash32(len, val, seed)

# 32-bit hash (calling C version in Base)
function h32c(s::MaybeSub, h::UInt) ; h += mhseed32 ; mh32c(s, h) + h ; end

# 32-bit hash (calling Julia version)
function h32j(s::MaybeSub, h::UInt) ; h += mhseed32 ; mh32j(s, h) + h ; end

mh128c(siz::Int, ptr::Ptr{UInt8}, seed=0) =
    ccall(:memhash_seed, UInt64, (Ptr{UInt8}, Csize_t, UInt32), ptr, siz, seed)
mh128c(s::MaybeSub, h=0) = mh128c(sizeof(s), pointer(s), h%UInt32)

mh128j4(s::MaybeSub, h=0) = mmhash128_4(sizeof(s), pointer(s), h%UInt32)
mh128j8a(s::MaybeSub, h=0) = mmhash128_8_a(sizeof(s), pointer(s), h%UInt32)
mh128j8u(s::MaybeSub, h=0) = mmhash128_8_u(sizeof(s), pointer(s), h%UInt32)
mh128j8c(s::MaybeSub, h=0) = mmhash128_8_c(s, h%UInt32)

# 128-bit MurmurHash (either 32-bit or 64-bit implementation, C version in Base), lower 64-bits
mh64c(siz::Int, ptr::Ptr{UInt8}, seed=0) =
    ccall(:memhash_seed, UInt, (Ptr{UInt8}, Csize_t, UInt32), ptr, siz, seed%UInt32)
mh64c(str::MaybeSub, seed=0) = mh64c(sizeof(str), pointer(str), seed%UInt32)

# 128-bit MurmurHash (Julia version)
mh128j(str::MaybeSub, seed=0) = mmhash128_a(sizeof(str), pointer(str), seed%UInt32)
mh128j(len::Int, val::Unsigned, seed=0) = mmhash128_a(len, val, seed%UInt32)
mh128j_c(str::AbstractString, seed=0) = mmhash128_c(str, seed%UInt32)

# 64-bit hash (calling 128-bit MurmurHash C version in Base)
function h64c(s::MaybeSub, h::UInt) ; h += mhseed64 ; mh64c(s, h%UInt32) + h ; end

# 64-bit hash (calling 128-bit Julia version)
function h64j(s, h::UInt) ; h += mhseed64 ; last(mh128j(s, h%UInt32)) + h ; end

load_u64(p) = unsafe_load(reinterpret(Ptr{UInt64}, pointer(p1)))

const sizp2 = sizeof(p2)
const unsp2 = load_u64(p2)

@testset "MurmurHash3" begin
    @testset "Aligned vs unaligned" begin
        @test mh128j_c(p1) == mh128j(p2)
    end
    @testset "32-bit MurmurHash" begin
        @test mh32j(p1) == mh32c(p1)
        @test mh32j(p2) == mh32c(p2)
        @test mh32j(sizp2, unsp2) == mh32c(p2)
    end
    @testset "128-bit MurmurHash" begin
        @test last(mh128j(sizp2, unsp2)) == mh64c(p1)
        @test last(mh128j(p1)) == mh64c(p1)
        @test last(mh128j_c(p2)) == mh64c(p2)
        @test last(mh128j_c(p2))  == mh64c(p1)
    end
end

print("32c:   "); dump(mh32c(p2))
print("32j:   "); dump(mh32j(p2))
print("32u:   "); dump(mh32j(sizp2, unsp2))
print("128c:  "); dump(mh128c(p2))
print("128j:  "); dump(mh128j(p2))
print("128u:  "); dump(mh128j(sizp2, unsp2))
print("128j4: "); dump(mh128j4(p2))
print("128ja: "); dump(mh128j8a(p2))
print("128ju: "); dump(mh128j8u(p2))
print("128jc: "); dump(mh128j8c(p2))
