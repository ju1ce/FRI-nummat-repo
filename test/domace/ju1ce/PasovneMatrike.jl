using Test

@testset "matrike" begin
    A = [1 2 3 ; 4 5 6 ; 7 8 9]
    Ap = PasovnaMatrika([1, 5, 9],[[4,8],[7]],[[2,6],[3]])
    
    eps = 1e-5
    
    @test norm(A-Ap)<eps
    
    A[1,1] = 10
    A[2,3] = 20
    A[3,1] = 30
    
    Ap[1,1] = 10
    Ap[2,3] = 20
    Ap[3,1] = 30
     
    @test norm(A-Ap)<eps
    
    A = [1 2 3 ; 0 5 6 ; 0 0 9]
    Ap = ZgornjePasovnaMatrika([1, 5, 9],[[2,6],[3]])
    
    @test norm(A-Ap)<eps
    
    A[1,1] = 10
    A[2,3] = 20
    
    Ap[1,1] = 10
    Ap[2,3] = 20
    
    @test norm(A-Ap)<eps
    
    A = [1 0 0 ; 4 5 0 ; 7 8 9]
    Ap = SpodnjePasovnaMatrika([1, 5, 9],[[4,8],[7]])
    
    eps = 1e-5
    
    @test norm(A-Ap)<eps
    
    A[1,1] = 10
    A[3,1] = 30
    
    Ap[1,1] = 10
    Ap[3,1] = 30
     
    @test norm(A-Ap)<eps
    
end

@testset "množenje" begin
    A = [1 2 3 ; 4 5 6 ; 7 8 9]
    Ap = PasovnaMatrika([1, 5, 9],[[4,8],[7]],[[2,6],[3]])
    
    eps = 1e-5
    
    v = [1,2,3]
    
    res = A*v
    resp = Ap*v
    
    @test norm(res-resp) < eps
    
    A = [1 2 3 ; 0 5 6 ; 0 0 9]
    Ap = ZgornjePasovnaMatrika([1, 5, 9],[[2,6],[3]])

    v = [1,2,3]
    
    res = A*v
    resp = Ap*v
    
    @test norm(res-resp) < eps
    
    A = [1 0 0 ; 4 5 0 ; 7 8 9]
    Ap = SpodnjePasovnaMatrika([1, 5, 9],[[4,8],[7]])
    
    v = [1,2,3]
    
    res = A*v
    resp = Ap*v
    
    @test norm(res-resp) < eps
    
    Ap = PasovnaMatrika(rand(10),[rand(9),rand(8),rand(7),rand(6)],[rand(9),rand(8),rand(7)])
    A = copy(Ap)
    
    v = rand(10)
    
    @test norm(Ap-A) < eps           #sanity check
    @test typeof(A) != typeof(Ap)
    
    res = A*v
    resp = Ap*v
    
    @test norm(res-resp) < eps
    
    Ap = ZgornjePasovnaMatrika(rand(10),[rand(9),rand(8),rand(7)])
    A = copy(Ap)
    
    v = rand(10)
    
    @test norm(Ap-A) < eps           #sanity check
    @test typeof(A) != typeof(Ap)
    
    res = A*v
    resp = Ap*v
    
    @test norm(res-resp) < eps
    
    Ap = SpodnjePasovnaMatrika(rand(10),[rand(9),rand(8),rand(7),rand(6)])
    A = copy(Ap)
    
    v = rand(10)
    
    @test norm(Ap-A) < eps           #sanity check
    @test typeof(A) != typeof(Ap)
    
    res = A*v
    resp = Ap*v
    
    @test norm(res-resp) < eps
    
end

@testset "deljenje" begin
    Ap = PasovnaMatrika(rand(10),[rand(9),rand(8),rand(7),rand(6)],[rand(9),rand(8),rand(7)])
    A = copy(Ap)
    v = rand(10)
    
    eps = 1e-5
    
    @test norm(Ap-A) < eps           #sanity check
    @test typeof(A) != typeof(Ap)
    
    res = A\v
    resp = Ap\v
    
    @test norm(res-resp) < eps
    
    Ap = ZgornjePasovnaMatrika(rand(10),[rand(9),rand(8),rand(7)])
    A = copy(Ap)
    v = rand(10)
    
    eps = 1e-5
    
    @test norm(Ap-A) < eps           #sanity check
    @test typeof(A) != typeof(Ap)
    
    res = A\v
    resp = Ap\v
    
    @test norm(res-resp) < eps
    
    Ap = SpodnjePasovnaMatrika(rand(10),[rand(9),rand(8),rand(7),rand(6)])
    A = copy(Ap)
    v = rand(10)

    eps = 1e-5
    
    @test norm(Ap-A) < eps           #sanity check
    @test typeof(A) != typeof(Ap)
    
    res = A\v
    resp = Ap\v
    
    @test norm(res-resp) < eps
    
end

@testset "lu" begin
    Ap = PasovnaMatrika(rand(10),[rand(9),rand(8),rand(7),rand(6)],[rand(9),rand(8),rand(7)])

    for i = 1:length(Ap.d)         #ensure diagonal dominance
         Ap[i,i] = sum(Ap[i,:])+1
    end

    A = copy(Ap)
    
    eps = 1e-5
    
    @test norm(Ap-A) < eps           #sanity check
    @test typeof(A) != typeof(Ap)
    
    Lp, Up = lu(Ap)
    L, U = lu(A)
    
    @test norm(L-Lp) < eps
    @test norm(U-Up) < eps
    
    @test norm(Lp*Up - Ap) < eps
    
end
