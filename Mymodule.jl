module Mymodule
export ∂x
export ∂y
export ∂z
export InterΣᶜᶜᶜ
export ΣΣ
export trΣ
export Σᵢ₃ᶜᶜᶜ 
######partial x
@inline function ∂x(u,Δx::Float64)
            Nx=size(u,1);Ny=size(u,2);Nz=size(u,3);Nt=size(u,4)
            ∂x_3d=zeros(Nx,Ny,Nz,Nt)
@inbounds    for k in 1:Nz,j in 1:Ny,i in 1:Nx
                        i₀=i+1
                        if i==Nx
                            i₀=1
                        end
                        ∂x_3d[i,j,k,:]=(u[i₀,j,k,:]-u[i,j,k,:])/Δx
            end
        return ∂x_3d
        end
######
##### partial y
@inline function ∂y(u,Δy::Float64)
            Nx=size(u,1);Ny=size(u,2);Nz=size(u,3);Nt=size(u,4)
            ∂y_3d=zeros(Nx,Ny,Nz,Nt)
@inbounds    for k in 1:Nz,j in 1:Ny,i in 1:Nx
                        j₀=j+1
                        if j==Ny
                            j₀=1
                        end
                        ∂y_3d[i,j,k,:]=(u[i,j₀,k,:]-u[i,j,k,:])/Δy
            end
        return ∂y_3d
        end
######
##### partial z
@inline function ∂z(u,Δz::Float64)
            Nx=size(u,1);Ny=size(u,2);Nz=size(u,3);Nt=size(u,4)
            ∂z_3d=zeros(Nx,Ny,Nz,Nt)
@inbounds    for k in 1:Nz-1,j in 1:Ny,i in 1:Nx
                        k₀=k+1
                        # if k==Nz
                        #     k₀=1
                        # end
                        ∂z_3d[i,j,k,:]=(u[i,j,k₀,:]-u[i,j,k,:])/Δz
            end
        return ∂z_3d
        end
#######
####Interpolation for Σ
@inline function InterΣᶜᶜᶜ(u,v,w,Δx,Δy,Δz)
#############对需要在z方向求偏导的参数做插值，避免底部固边界的插值问题
            uᵢ=zeros(size(w,1),size(w,2),size(w,3),size(w,4))
            vᵢ=zeros(size(w,1),size(w,2),size(w,3),size(w,4))
@inbounds    for k in 1:size(u,3),j in 1:size(u,2),i in 1:size(u,1)
                        if k != 1
                            uᵢ[i,j,k,:]=(u[i,j,k-1,:]+u[i,j,k,:])./2
                            vᵢ[i,j,k,:]=(v[i,j,k-1,:]+v[i,j,k,:])./2
                        elseif k==size(u,3)
                            uᵢ[:,:,k,:]=u[:,:,k,:]
                            vᵢ[:,:,k,:]=v[:,:,k,:]
                        end
            end
####################这里∂u∂z等函数替换成插值后的函数方便后续处理                
            ∂u∂y=∂y(u,Δy);∂u∂z=∂z(uᵢ,Δz)
            ∂v∂x=∂x(v,Δx);∂v∂z=∂z(vᵢ,Δz)
            ∂w∂x=∂x(w,Δx);∂w∂y=∂y(w,Δy)
            Σ=zeros(size(u,1),size(u,2),size(u,3),size(u,4))
@inbounds    for k in 1:size(u,3),j in 1:size(u,2),i in 1:size(u,1)
#####################应对周期性边界条件
                        i₊₁=i+1; i₋₁=i-1; j₊₁=j+1; j₋₁=j-1; k₊₁=k+1; k₋₁=k-1
                        if i₊₁==size(u,1)+1
                            i₊₁=1
                        end
                        if i₋₁==0
                            i₋₁=size(u,1)
                        end
                        if j₊₁==size(u,2)+1
                            j₊₁=1
                        end
                        if j₋₁==0
                            j₋₁=size(u,2)
                        end
######################求Σ
                        Σ[i,j,k,:]=((∂u∂y[i,j,k,:]+∂u∂y[i₊₁,j,k,:]+∂u∂y[i,j₋₁,k,:]+∂u∂y[i₊₁,j₋₁,k,:])/4+(∂v∂x[i,j,k,:]+∂v∂x[i₋₁,j,k,:]+∂v∂x[i,j₊₁,k,:]+∂v∂x[i₋₁,j₊₁,k,:])/4
                        +(∂u∂z[i,j,k,:]+∂u∂z[i₊₁,j,k,:])/2+(∂w∂x[i,j,k,:]+∂w∂x[i₋₁,j,k,:]+∂w∂x[i,j,k₊₁,:]+∂w∂x[i₋₁,j,k₊₁,:])/4
                        +(∂v∂z[i,j,k,:]+∂v∂z[i,j₊₁,k,:])/2+(∂w∂y[i,j,k,:]+∂w∂y[i,j₋₁,k,:]+∂w∂y[i,j,k₊₁,:]+∂w∂y[i,j₋₁,k₊₁,:])/4)
            end
        return Σ
    end
##################Σ²
@inline function ΣΣ(u,v,w,Δx,Δy,Δz)
                return InterΣᶜᶜᶜ(u,v,w,Δx,Δy,Δz).^2
            end


@inline function trΣ(u,v,w,Δx,Δy,Δz)
                trΣ=zeros(size(u,1),size(u,2),size(u,3),size(u,4))
                ∂u∂x=∂x(u,Δx);∂v∂y=∂y(v,Δy);∂w∂z=∂z(w,Δz);
    @inbounds   for k in 1:size(u,3)-1,j in 1:size(u,2),i in 1:size(u,1)
                    trΣ[i,j,k,:]=∂u∂x[i,j,k,:]+∂v∂y[i,j,k,:]+∂w∂z[i,j,k,:]
                end
        return trΣ
    end
#################Σᵢ₃ for Shear production
@inline function Σᵢ₃ᶜᶜᶜ(u,v,w,Δx,Δy,Δz)
#############对需要在z方向求偏导的参数做插值，避免底部固边界的插值问题
            uᵢ=zeros(size(w,1),size(w,2),size(w,3),size(w,4))
            vᵢ=zeros(size(w,1),size(w,2),size(w,3),size(w,4))
@inbounds    for k in 2:size(u,3),j in 1:size(u,2),i in 1:size(u,1)
                    if k != 1
                        uᵢ[i,j,k,:]=(u[i,j,k-1,:]+u[i,j,k,:])./2
                        vᵢ[i,j,k,:]=(v[i,j,k-1,:]+v[i,j,k,:])./2
                    elseif k==size(u,3)
                        uᵢ[:,:,k,:]=u[:,:,k,:]
                        vᵢ[:,:,k,:]=v[:,:,k,:]
                    end
            end
####################这里∂u∂z等函数替换成插值后的函数方便后续处理                
            ∂u∂z=∂z(uᵢ,Δz);∂v∂z=∂z(vᵢ,Δz)
            ∂w∂x=∂x(w,Δx);∂w∂y=∂y(w,Δy)
            Σᵢ₃=zeros(size(u,1),size(u,2),size(u,3),size(u,4))
 @inbounds  for k in 1:size(u,3),j in 1:size(u,2),i in 1:size(u,1)
#####################应对周期性边界条件
                    i₊₁=i+1; i₋₁=i-1; j₊₁=j+1; j₋₁=j-1; k₊₁=k+1; k₋₁=k-1
                    if i₊₁==size(u,1)+1
                        i₊₁=1
                    end
                    if i₋₁==0
                        i₋₁=size(u,1)
                    end
                    if j₊₁==size(u,2)+1
                        j₊₁=1
                    end
                    if j₋₁==0
                        j₋₁=size(u,2)
                    end
######################求Σᵢ₃
                    Σᵢ₃[i,j,k,:]=((∂u∂z[i,j,k,:]+∂u∂z[i₊₁,j,k,:])/2+(∂w∂x[i,j,k,:]+∂w∂x[i₋₁,j,k,:]+∂w∂x[i,j,k₊₁,:]+∂w∂x[i₋₁,j,k₊₁,:])/4
                                +(∂v∂z[i,j,k,:]+∂v∂z[i,j₊₁,k,:])/2+(∂w∂y[i,j,k,:]+∂w∂y[i,j₋₁,k,:]+∂w∂y[i,j,k₊₁,:]+∂w∂y[i,j₋₁,k₊₁,:])/4)
            end
        return Σᵢ₃
    end
end #end for module
