function plot_3D(x)
    n = size(x)
    let
        vol = [x[ix,iy,iz] for ix in 1:n[1], iy in 1:n[2], iz in 1:n[3]]
        fig, ax, _ = volume(1:n[1], 1:n[2], 1:n[3], vol, colormap = :plasma, colorrange = (minimum(vol), maximum(vol)),
            figure = (; resolution = (800,800)),  
            axis=(; type=Axis3, perspectiveness = 0.5,  azimuth = 7.19, elevation = 0.57,  
                aspect = (1,1,1)))
    
        fig
    end
end