% Post-processing MATLAB script for extracting 1-D properties
% at certain horizontal points
% Modified from Nakagawa scripts
%

clear

directory = '../NCFMAS-R/';
file_name = 'perplex_binary_2d1e8_R104';
horizontal_points=[1 128 256 384 512 640 758 896];


file_stem = [directory file_name];

for frame = 37:1:37
% Frame number: Please refer the final two digit in your data file
frame

if frame<10000
  framestring = num2str(10000+frame);
  framestring(1) = '0';
end

% load spatial fields: 4D arrays (x,y,z,b)

[theta, phi, z, vs, time] = ReadStag3Dpjt(directory, file_name, frame, 's-velocity');
[theta, phi, z, vb, time] = ReadStag3Dpjt(directory, file_name, frame, 'b-velocity');
[theta, phi, z, ps, time] = ReadStag3Dpjt(directory, file_name, frame, 'p-velocity');
[theta, phi, z,  rho, time] = ReadStag3Dpjt(directory, file_name, frame, 'density');
[theta, phi, z,  T, time] = ReadStag3Dpjt(directory, file_name, frame, 'temperature');
[theta, phi, z,  c, time] = ReadStag3Dpjt(directory, file_name, frame, 'composition');

[nt, np, nz]=size(vs);
nz

numel(horizontal_points)
vs_arr=zeros(numel(horizontal_points), nz);
vb_arr=zeros(numel(horizontal_points), nz);
vp_arr=zeros(numel(horizontal_points), nz);
rho_arr=zeros(numel(horizontal_points), nz);
T_arr=zeros(numel(horizontal_points), nz);
c_arr=zeros(numel(horizontal_points), nz);

for k = 1:(numel(horizontal_points)-1)
  p = horizontal_points(k)
  vs_arr(k,:) = squeeze(vs(1,p,:));
  vb_arr(k,:) = squeeze(vb(1,p,:));
  vp_arr(k,:) = squeeze(ps(1,p,:)); % sqrt(vb_arr(k).^2+4.0/3.0*vs_arr(k).^2);

  rho_arr(k,:) = squeeze(rho(1,p,:));
  T_arr(k,:) = squeeze(T(1,p,:));
  c_arr(k,:) = squeeze(c(1,p,:));
end

for k = 1:numel(horizontal_points)-1
  p = horizontal_points(k)

  tname=strcat(file_stem,'_all_',sprintf('%05d',p),'_', framestring,'.dat');
  fid=fopen(tname,'w');
  for l = 1:nz
     fprintf(fid,'%f %f %f %f %f %f %f\n', (1.0-z(l))*2890, T_arr(k,l)*2500, c_arr(k,l), rho_arr(k,l), vs_arr(k,l), vb_arr(k,l), vp_arr(k,l));
  end
end

end
