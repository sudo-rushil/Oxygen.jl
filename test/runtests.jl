using Oxygen, Test

@testset "foo" begin
    x, y = 5, 7
    @test foo(x, y) == 7
    x = "blah"
    @test_throws MethodError foo(x, y)
end

@testset "bar" begin
    z = 4.
    @test bar(z) == 1.
end

@testset "baz" begin
    n = 5
    @test baz(n) == 5
    n = "blah"
    @test_throws MethodError baz(n)
end
