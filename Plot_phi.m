clear all;
lighting gouraud;
%vidObj = VideoWriter('movie.avi');
%vidObj.Quality = 100;
%vidObj.FrameRate = 10;
%open(vidObj);

%Loading the file[=
%for index=0:2000:200000
    clear figure
    %camlight off
%    frame=1;
    index = 2000;
    filename = sprintf('%s%d','prof_',index);
    %disp(filename);
    data = load(filename);
    nx   = 128;
    ny   = 128;
    nz   = 128;
    phi  = rand(nx,ny,nz);

    for i = 0:127
        for j = 0:127
            for k = 0:127
                l = (k+nz*(j+i*ny))+1;
                phi(i+1,j+1,k+1) = data(l,4);
                %disp(l);
            end
        end
    end

    data = smooth3(phi,'box',5);
    [x,y,z] = meshgrid(1:nx,1:ny,1:nz);
    axis([1 nx 1 ny 1 nz]);
    p1 = patch(isosurface(x,y,z,phi,0.5));
    isonormals(phi,p1);
    p1.FaceColor = 'yellow';
    p1.EdgeColor = 'none';
    
    p2 = patch(isocaps(phi,0.5));
    p2.FaceColor = 'yellow';
    p2.EdgeColor = 'none';
    
    p3 = patch(isosurface(x,y,z,1.0-phi,0.5));
    isonormals(phi,p1);
    p3.FaceColor = [1.0 0.27 0.0];
    p3.EdgeColor = 'none';
    
    p4 = patch(isocaps(1.0-phi,0.5));
    p4.FaceColor = [1.0 0.27 0.0];
    p4.EdgeColor = 'none';
    
    %p3 = patch(isocaps(phi, -0.1));
    %p3.FaceColor = 'interp';
    %p3.EdgeColor = 'none';
    %colormap rgb;
    %isosurface(x,y,z,phi,0.5);
    %lightangle(-45,30)
    %p.FaceLighting = 'gouraud';
    %p.AmbientStrength = 0.3;
    %p.DiffuseStrength = 0.8;
    %p.SpecularStrength = 0.9;
    %p.SpecularExponent = 25;
    %p.BackFaceLighting = 'unlit';
    daspect([1 1 1]);
    view(-150,30);
    camlight left;
    camlight right;
    lighting gouraud;
    axis square;
    %writeVideo(vidObj,getframe(gcf));
    %delete(findall(gcf,'Type','light'));
%end
%close(vidObj);
