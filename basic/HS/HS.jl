using LinearAlgebra
using Plots

#---------------------------------------------------------------------
a=1
α=acos(exp(-a))
C=1/2
x=collect(-2:0.01:2)
dy=zeros(length(x))
for i in eachindex(x)

    dy[i]=exp(-a*x[i]^2)-C*exp(1im*α*x[i])-C*exp(-1im*α*x[i])

end
scatter(x,dy,label="Two component ")

savefig("E:\\桌面\\JuliaDQMC\\basic\\HS\\1.png")

#---------------------------------------------------------------------

a=1
α=sqrt(a)
γ=[1+sqrt(6)/3 1+sqrt(6)/3 1-sqrt(6)/3 1-sqrt(6)/3]
η=[sqrt(2*(3-sqrt(6))) -sqrt(2*(3-sqrt(6))) sqrt(2*(3+sqrt(6))) -sqrt(2*(3+sqrt(6)))]
C=1/4

x=collect(-1:0.01:1)
dy=zeros(length(x))*1im

for i in eachindex(x)
    dy[i]=exp(-a*x[i]^2)
    for j in 1:4
        dy[i]-=C*γ[j]*exp(1im*sqrt(a)*η[j]*x[i] )
    end
end

scatter(x,abs.(dy),label="Four component ")
savefig("E:\\桌面\\JuliaDQMC\\basic\\HS\\2.png")

#---------------------------------------------------------------------

# 四分量在O=±1随a的误差
a=collect(0.0001:0.001:1)
y=2*exp.(-a)-(1+sqrt(6)/3)*cos.(sqrt.( 2*(3-sqrt(6))*a )) - (1-sqrt(6)/3)*cos.(sqrt.( 2*(3+sqrt(6))*a ))
plot(a,y)

# attactive-U 四分量η确定
# x=collect(0:0.01:2)
# y=collect(0:0.01:2)
# a=0.01




# 1+sqrt(6)/3
# 1-sqrt(6)/3


#---------------------------------------------------------------------

a=1
α=sqrt(a)
γ=[1+sqrt(6)/3 1+sqrt(6)/3 1-sqrt(6)/3 1-sqrt(6)/3]
η=[sqrt(2*(3-sqrt(6))) -sqrt(2*(3-sqrt(6))) sqrt(2*(3+sqrt(6))) -sqrt(2*(3+sqrt(6)))]
C=1/4

x=collect(-1:0.01:1)
dy=zeros(length(x))*1im

for i in eachindex(x)
    dy[i]=exp(a*x[i]^2)
    for j in 1:4
        dy[i]-=C*γ[j]*exp(-sqrt(a)*η[j]*x[i] )
    end
end

scatter(x,abs.(dy),label="Four component ")
savefig("E:\\桌面\\JuliaDQMC\\basic\\HS\\3.png")


