
Random.seed!(25)

include("./helpers/utils.jl")

# set up.
fs = 14005.602240896402
SW = 20.0041938620844
ν_0ppm = 10656.011933076665
hz2ppmfunc = uu->(uu - ν_0ppm)*SW/fs
ppm2hzfunc = pp->(ν_0ppm + pp*fs/SW)

u_rad_test = ppm2hzfunc(3.3984) *2*π

L = 20
ds = randn(L)
ξs = randn(L)

real_setp, imag_setp, _ = loadalphaglucose()

# timing.
out2r = real_setp(u_rad_test - ds[l], ξs[l])
out2i = imag_setp(u_rad_test - ds[l], ξs[l])

# the following are 50 ns each.
@btime real_setp($ds[l], $ξs[l]) setup=(l=rand(1:length($ds))) # 50.
@btime real_setp(x,y) setup=(x=rand(),y=randn()) # 62ns.

y = rand()
@btime real_setp(x,$y) setup=(x=rand()) # 62ns.

c1 = rand()
c2 = rand()
@btime real_setp($c1,$c2)


# subtraction at input.
@btime x+1 setup=(x=randn())
@btime real_setp($u_rad_test - $ds[l], $ξs[l]) setup=(l=rand(1:length($ds)))


@assert 1==2

ds = FM.ds
@btime u_rad_test - ds[3]


x = randn()
y = randn()
z = u_rad_test - ds[3]
@btime real_setp(z, x)
@btime x - y 



@btime x-y setup=(x=rand(),y=randn()) # 17ns
@btime real_setp(x,y) setup=(x=ppm2hzfunc(rand()+3)*2*π,y=randn()) # 61 ns
@btime real_setp(x,y) setup=(x=rand(),y=randn()) # 62ns.
@btime real_setp(u_rad_test - FM.ds[l], FM.ξs[l]) setup=(l=rand(1:length(FM.ds))) # 227ns.

@btime real_setp(FM.ds[l], FM.ξs[l]) setup=(l=rand(1:length(FM.ds))) # 194.

# idea: f.(flat) to speed up.