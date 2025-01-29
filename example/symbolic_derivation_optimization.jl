using Symbolics
using LinearAlgebra
@variables ori_x ori_y dire_x dire_y  theta hh2 bb  t1 ang1 dist
@variables AB[2] dist2 Tₑ² Tₑ[2] pp[2]  dot_AB_P dot_AB_Tₑ dot_AB_D dot_D_Tₑ dot_D_P
@variables N1 N2
@variables ∂²f∂θ∂t ∂²f∂θ² ∂²f∂t²
@variables ∂²f∂s₁²  ∂²f∂s₂²  ∂²f∂s₃²
@variables ∂²f∂s₁∂t  ∂²f∂s₂∂θ  ∂²f∂s₃∂θ
@variables s1 s2 s3 lambda1 lambda2 lambda3 p_barrier upper lower
@variables ∂f∂t ∂f∂θ ∂f∂s1 ∂f∂s2 ∂f∂s3
@variables C1
ori=[ori_x;ori_y]
dire=[dire_x;dire_y]
p=[cos(theta);bb*sin(theta)]

V=ori+dire*t1-p

fun=V'*V-hh2

fun2=fun*fun

fun2

Dt1=Differential(t1)
Dtheta=Differential(theta)
fun2
_D_t=expand_derivatives(Dt1(fun2))

D_t= _D_t |> x-> substitute(x, [ori_x=> AB[1]-(dire_x*t1-cos(theta)), ori_y=> AB[2]-(dire_y*t1-bb*sin(theta)),-hh2=>dist2-AB[1]^2-AB[2]^2]) |> simplify
_D_theta=expand_derivatives(Dtheta(fun2))
D_theta=_D_theta |> x-> substitute(x, [ori_x=> AB[1]-(dire_x*t1-cos(theta)), ori_y=> AB[2]-(dire_y*t1-bb*sin(theta)),-hh2=>dist2-AB[1]^2-AB[2]^2]) |> simplify


D_t2=expand_derivatives(Dt1(_D_t))  |> x-> substitute(x, [ori_x=> AB[1]-(dire_x*t1-cos(theta)), ori_y=> AB[2]-(dire_y*t1-bb*sin(theta)),-hh2=>dist2-AB[1]^2-AB[2]^2,dire_x^2=>1-dire_y^2]) |> simplify |>
x-> substitute(x, (2AB[1]*dire_x + 2AB[2]*dire_y) => 2*dot_AB_D) |> simplify
_D_theta2=expand_derivatives(Dtheta(_D_theta))  |> x-> substitute(x, [ori_x=> AB[1]-(dire_x*t1-cos(theta)), ori_y=> AB[2]-(dire_y*t1-bb*sin(theta)),-hh2=>dist2-AB[1]^2-AB[2]^2]) |> simplify

_coeff_dist2=Symbolics.coeff(_D_theta2,dist2)

coeff_dist2 = 2(2AB[1]*cos(theta) + 2(sin(theta)^2) + 2AB[2]*bb*sin(theta) + 2(bb^2)*(cos(theta)^2)) |> x-> substitute(x, [sin(theta)^2=>Tₑ²-bb^2*cos(theta)^2]) |> expand |> simplify |>
x-> substitute(x, 4*AB[1]*cos(theta)=> 4*dot_AB_P-4*AB[2]*bb*sin(theta)) |> expand |> simplify

D_theta2= _D_theta2 |> x-> substitute(x, [2(2AB[1]*cos(theta) + 2(sin(theta)^2) + 2AB[2]*bb*sin(theta) + 2(bb^2)*(cos(theta)^2))*dist2=>coeff_dist2*dist2,(2AB[1]*sin(theta) - 2AB[2]*bb*cos(theta))=>2*dot_AB_Tₑ]) |> simplify


_D_t_theta=expand_derivatives(Dtheta(_D_t)) |> x-> substitute(x, [ori_x=> AB[1]-(dire_x*t1-cos(theta)), ori_y=> AB[2]-(dire_y*t1-bb*sin(theta)),-hh2=>dist2-AB[1]^2-AB[2]^2,dire_x^2=>1-dire_y^2]) |> simplify |>
x-> substitute(x, (2AB[1]*dire_x + 2AB[2]*dire_y) => 2*dot_AB_D) |> simplify
_D_theta_t=expand_derivatives(Dt1(_D_theta))

D_t_theta= _D_t_theta |> x-> substitute(x, [(2AB[1]*sin(theta) - 2AB[2]*bb*cos(theta))=> 2*dot_AB_Tₑ]) |> simplify |>
x-> substitute(x,2*(2dire_x*sin(theta) - 2bb*dire_y*cos(theta))*dist2=>4*dot_D_Tₑ*dist2) |> simplify

@info "D_t2" D_t2
@info "D_theta2" D_theta2
@info "D_t_theta" D_t_theta

@info "D_theta" D_theta
@info "D_t" D_t

D_theta2*D_t2  |>expand
-D_t_theta^2 |> expand
denn=D_t2*D_theta2-D_t_theta*D_t_theta |> expand #|> x-> x/dist2 |> expand |> simplify
dist_denn_coeff=Symbolics.coeff(denn,dist2)

denn0=denn-dist_denn_coeff*dist2 |> expand |> simplify

-denn0/dist_denn_coeff |> expand |> simplify

Tange=[sin(theta);-bb*cos(theta)]

ABTe=dot(V,Tange) |> expand |> simplify
ABP=dot(V,[cos(theta);bb*sin(theta)]) |> expand |> simplify
ABTe2=ABTe^2 |> expand |> simplify

_h_theta=dot_AB_Tₑ*Dtheta(ABTe)*2 |> expand_derivatives |>
x-> substitute(x, [ori_x*cos(theta)=> dot_AB_P-(dire_x*t1*cos(theta)+dire_y*bb*t1*sin(theta)+ori_y*bb*sin(theta))]) |> simplify |>
x-> substitute(x,(-(1//1) + bb^2)=> -C1) |> simplify
_h_t1=2*dot_AB_Tₑ*Dt1(ABTe) |> expand_derivatives |> simplify |>
x-> substitute(x, (dire_x*sin(theta) - bb*dire_y*cos(theta))=> dot_D_Tₑ) |> simplify


DTₑ=dot([dire_x;dire_y],Tange) |> expand |> simplify

__a1=2*Dtheta(DTₑ)*dot_AB_Tₑ+2*dot_D_Tₑ*Dtheta(ABTe) |> expand_derivatives |> simplify |>
x-> substitute(x, [ori_x*cos(theta)=> dot_AB_P-(dire_x*t1*cos(theta)+dire_y*bb*t1*sin(theta)+ori_y*bb*sin(theta))]) |>
x-> substitute(x, (dire_x*cos(theta) + bb*dire_y*sin(theta))=>dot_D_P) |> simplify |>
x-> substitute(x, (dot_AB_P + (-(1//1) + bb^2)*cos(2theta))=>(dot_AB_P-C1*cos(2theta))) |> simplify

__a2=2*Dt1(ABP)*dot_AB_Tₑ+2*(dot_AB_P-C1*cos(2theta))*Dt1(ABTe) |> expand_derivatives |> simplify |>
x-> substitute(x, [ori_x*cos(theta)=> dot_AB_P-(dire_x*t1*cos(theta)+dire_y*bb*t1*sin(theta)+ori_y*bb*sin(theta))]) |>
x-> substitute(x, (dire_x*sin(theta) - bb*dire_y*cos(theta))=>dot_D_Tₑ) |> simplify |>
x-> substitute(x, (dire_x*cos(theta) + bb*dire_y*sin(theta))=>dot_D_P) |> simplify



__a1-__a2


2dot_AB_Tₑ*dot_D_P - 2(-1 + dot_AB_P + bb^2)*dot_D_Tₑ + 2(dot_AB_P + (-1+ bb^2)*cos(2theta))*dot_D_Tₑ - 2(dire_x*cos(theta) + bb*dire_y*sin(theta))*dot_AB_Tₑ

_h_theta2=2*Dtheta(ABP)*dot_AB_Tₑ+2*(dot_AB_P-C1*cos(2theta))*Dtheta(ABTe) |> expand_derivatives |>
x-> substitute(x, [ori_x*cos(theta)=> dot_AB_P-(dire_x*t1*cos(theta)+dire_y*bb*t1*sin(theta)+ori_y*bb*sin(theta))]) |>
x-> substitute(x, [- ori_x*sin(theta)=> -dot_AB_Tₑ+(dire_x*t1*sin(theta)-dire_y*bb*t1*cos(theta)-ori_y*bb*cos(theta))]) |> simplify |>
x-> substitute(x, (-dot_AB_Tₑ + sin(2theta) - (bb^2)*sin(2theta))=> -(dot_AB_Tₑ-C1*sin(2theta))) |>
x-> substitute(x, (dot_AB_P + (-(1//1) + bb^2)*cos(2theta))=>(dot_AB_P-C1*cos(2theta))) |> simplify

_h_t12=2*(Dt1(ABTe)*dot_D_Tₑ+Dt1(DTₑ)*dot_AB_Tₑ) |> expand_derivatives |> simplify |>
x-> substitute(x, (dire_x*sin(theta) - bb*dire_y*cos(theta))=>dot_D_Tₑ) |> simplify


@info "∂f∂t: $_h_t1"
@info "∂f∂θ: $_h_theta"
@info "∂f∂t∂θ: $__a1"
@info "∂f∂θ∂t: $__a2"
@info "∂f∂θ2: $_h_theta2"
@info "∂f∂t2: $_h_t12"




ss1=s1#abs(s1)
ss2=s2#abs(s2)
ss3=s3#abs(s3)

h1=(t1-ss1)
g1=(theta-lower-ss2)
g2=(upper-theta-ss3)

h1²=h1^2
g1²=g1^2
g2²=g2^2
Ds1=Differential(s1)
Ds2=Differential(s2)
Ds3=Differential(s3)
_penalty=lambda1*h1+
lambda2*g1+
lambda3*g2+
p_barrier*(h1²+g1²+g2²)/2 |> expand

coeff_penalty1=Symbolics.coeff(_penalty,lambda1)
coeff_penalty2=Symbolics.coeff(_penalty,lambda2)
coeff_penalty3=Symbolics.coeff(_penalty,lambda3)
coeff_penalty_barrier=Symbolics.coeff(_penalty,p_barrier)  |>  simplify

vcat(coeff_penalty1,coeff_penalty2,coeff_penalty3,coeff_penalty_barrier)

∂²f1∂t²=Dt1(_penalty) |> expand_derivatives |> Dt1 |> expand_derivatives |> simplify
∂²f1∂θ²=Dtheta(_penalty) |> expand_derivatives |> Dtheta |> expand_derivatives |> simplify
∂²f1∂s1²=Ds1(_penalty) |> expand_derivatives |> Ds1 |> expand_derivatives |> simplify
∂²f1∂s2²=Ds2(_penalty) |> expand_derivatives |> Ds2 |> expand_derivatives |> simplify
∂²f1∂s3²=Ds3(_penalty) |> expand_derivatives |> Ds3 |> expand_derivatives |> simplify
∂²f1∂s1∂t=Ds1(Dt1(_penalty)) |> expand_derivatives |> simplify
∂²f1∂s2∂θ=Ds2(Dtheta(_penalty)) |> expand_derivatives |> simplify
∂²f1∂s3∂θ=Ds3(Dtheta(_penalty)) |> expand_derivatives |> simplify
∂f1∂θ=Dtheta(_penalty) |> expand_derivatives |> simplify
∂f1∂t=Dt1(_penalty) |> expand_derivatives |> simplify
∂f1∂s1=Ds1(_penalty) |> expand_derivatives |> simplify
∂f1∂s2=Ds2(_penalty) |> expand_derivatives |> simplify
∂f1∂s3=Ds3(_penalty) |> expand_derivatives |> simplify


DF=[∂f∂t+∂f1∂t; ∂f∂θ+∂f1∂θ ; ∂f1∂s1; ∂f1∂s2; ∂f1∂s3]

F=[∂²f∂t²+∂²f1∂t²      ∂²f∂θ∂t      ∂²f1∂s1∂t        0        0;
  ∂²f∂θ∂t       ∂²f∂θ²+∂²f1∂θ²        0         ∂²f1∂s2∂θ   ∂²f1∂s3∂θ;
  ∂²f1∂s1∂t       0            ∂²f1∂s1²        0        0;
  0              ∂²f1∂s2∂θ      0            ∂²f1∂s2²   0;
  0              ∂²f1∂s3∂θ      0            0        ∂²f1∂s3²]

invF=inv(F) |> expand |> expand

denom=det(F) |> expand |>
x-> (0,Symbolics.coeff(x, p_barrier^3),substitute(x,p_barrier^3=>0),0,0) |>
expand |> simplify|> x-> vcat(x...)


H=invF|> x-> substitute.(x, (-8(p_barrier^3)*(∂²f∂θ∂t^2) + (2p_barrier + ∂²f∂t²)*(-32(p_barrier^4) + 8(4p_barrier + ∂²f∂θ²)*(p_barrier^3)) - 2p_barrier*(-32(p_barrier^4) - 8(-4p_barrier - ∂²f∂θ²)*(p_barrier^3)))=>8(p_barrier^3)*∂²f∂t²*∂²f∂θ² - 8(p_barrier^3)*(∂²f∂θ∂t^2)) |> x-> @. expand(x) |> x-> @. simplify(x)

H1=invF|> x-> substitute.(x,  (-8(p_barrier^3)*(∂²f∂θ∂t^2) + (2p_barrier + ∂²f∂t²)*(-32(p_barrier^4) + 8(4p_barrier + ∂²f∂θ²)*(p_barrier^3)) - 2p_barrier*(-32(p_barrier^4) - 8(-4p_barrier - ∂²f∂θ²)*(p_barrier^3)))=>1) |> x-> @. expand(x) |> x-> @. simplify(x)


Ciff=H*DF |> x-> @. expand(x) |> x-> @. simplify(x) |>
x-> @. substitute(x,(∂²f∂t²*∂²f∂θ² - (∂²f∂θ∂t^2))=>C1).*p_barrier

Ciff[5]

Ciff |> x-> @. substitute(x,s1=>0) |> x-> @. substitute(x,s2=>0) |> x-> @. substitute(x,s3=>0)

Coeff4=Ciff |> x-> @. Symbolics.coeff(x,p_barrier^4)
Coeff3=Ciff |> x-> @. Symbolics.coeff(x,p_barrier^3)
Coeff2=Ciff |> x-> @. Symbolics.coeff(x,p_barrier^2)
Coeff0=Ciff |> x-> @. Symbolics.coeff(x,p_barrier^1)  |> x-> @. expand(x) |> x-> @. simplify(x)
Coeffn1=Ciff-Coeff0*p_barrier |> x-> @. expand(x) |> x-> @. simplify(x)

Coeff0s1 = @. Symbolics.coeff(Coeff0,s1) |> x-> @. substitute(x,(∂²f∂t²*∂²f∂θ² - (∂²f∂θ∂t^2))=>C1)
Coeff0s2 = @. Symbolics.coeff(Coeff0,s2) |> x-> @. substitute(x,(∂²f∂t²*∂²f∂θ² - (∂²f∂θ∂t^2))=>C1)
Coeff0s3 = @. Symbolics.coeff(Coeff0,s3) |> x-> @. substitute(x,(∂²f∂t²*∂²f∂θ² - (∂²f∂θ∂t^2))=>C1)
Coeff0t1=  @. Symbolics.coeff(Coeff0,t1) |> x-> @. substitute(x,(-∂²f∂t²*∂²f∂θ² + (∂²f∂θ∂t^2))=>-C1)


Coeff0theta= @. Symbolics.coeff(Coeff0,theta) |> x-> @. substitute(x,(-∂²f∂t²*∂²f∂θ² +(∂²f∂θ∂t^2))=>-C1) |>
x-> @. substitute(x,(∂²f∂t²*∂²f∂θ² - (∂²f∂θ∂t^2))=>C1)

Coeff0lower= @. Symbolics.coeff(Coeff0,lower) |> x-> @. substitute(x,(-∂²f∂t²*∂²f∂θ² + (∂²f∂θ∂t^2))=>-C1) |>
x-> @. substitute(x,(∂²f∂t²*∂²f∂θ² - (∂²f∂θ∂t^2))=>C1)

Coeffupper= @. Symbolics.coeff(Coeff0,upper) |> x-> @. substitute(x,(-∂²f∂t²*∂²f∂θ² + (∂²f∂θ∂t^2))=>-C1) |>
x-> @. substitute(x,(∂²f∂t²*∂²f∂θ² - (∂²f∂θ∂t^2))=>C1)

Coeff0nots= @. substitute(Coeff0, s1=>0) |> x-> @. substitute(x, s2=>0) |> x-> @. substitute(x, s3=>0) |>
x-> @. substitute(x, t1=>0) |> x-> @. substitute(x, theta=>0) |> x-> @. substitute(x, lower=>0) |> x-> @. substitute(x, upper=>0) |>
x-> @. substitute(x, ∂f∂t*∂²f∂θ²-∂f∂θ*∂²f∂θ∂t =>N1) |>
x-> @. substitute(x, ∂f∂θ*∂²f∂t²-∂f∂t*∂²f∂θ∂t =>N2) |>
x-> @. substitute(x,(∂f∂t*∂²f∂θ∂t - ∂f∂θ*∂²f∂t²)=>-N2)

Coeffn1_l1= @. Symbolics.coeff(Coeffn1,lambda1) |> x-> @. substitute(x,(∂²f∂t²*∂²f∂θ² - (∂²f∂θ∂t^2))=>C1) |>
x-> @. substitute(x,-(∂²f∂t²*∂²f∂θ² - (∂²f∂θ∂t^2))=>-C1)
Coeffn1_l2= @. Symbolics.coeff(Coeffn1,lambda2) |> x-> @. substitute(x,(∂²f∂t²*∂²f∂θ² - (∂²f∂θ∂t^2))=>C1) |>
x-> @. substitute(x,-(∂²f∂t²*∂²f∂θ² - (∂²f∂θ∂t^2))=>-C1)
Coeffn1_l3= @. Symbolics.coeff(Coeffn1,lambda3)  |> x-> @. substitute(x,(∂²f∂t²*∂²f∂θ² - (∂²f∂θ∂t^2))=>C1) |>
x-> @. substitute(x,-(∂²f∂t²*∂²f∂θ² - (∂²f∂θ∂t^2))=>-C1)

Coeffn1_l1+Coeffn1_l2+Coeffn1_l3

Coeff0s1.*s1+Coeff0s2.*s2+Coeff0s3.*s3+Coeff0t1.*t1+Coeff0theta.*theta+Coeff0lower.*lower+Coeffupper.*upper+
Coeff0nots+Coeffn1_l1.*lambda1/p_barrier+Coeffn1_l2.*lambda2/p_barrier+Coeffn1_l3.*lambda3/p_barrier


# C1 -> (∂²f∂t²*∂²f∂θ² - (∂²f∂θ∂t^2))
# N1 -> ∂f∂t*∂²f∂θ²-∂f∂θ*∂²f∂θ∂t
# N2 -> ∂f∂θ*∂²f∂t²-∂f∂t*∂²f∂θ∂t
