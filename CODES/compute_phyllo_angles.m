function [phi_all, theta_all, Kpos] = compute_phyllo_angles(Data,ssi)

nr_rad_samples = size(Data,1);
 nr_rad_lines = size(Data,2);
 nr_interleaves = size(Data,3);
delta_kr = 1/size(Data,1);
rad_pos(:,1) = -0.5:delta_kr:0.5-delta_kr;
                    
for interleaf = 0:nr_interleaves-1
                        for py_number = 0:nr_rad_lines-1
                            
% % % %                         %% koosh ----------------
%                             z = (py_number + 0.5 - nr_rad_lines) / nr_rad_lines;
%                             phi = sqrt(2* nr_rad_lines * pi / nr_interleaves)*asin(z) + interleaf * 2 * pi / nr_interleaves;

%                         %% Phyllotaxis -------------------------------
                            N = nr_interleaves * nr_rad_lines;%                             
                            n = py_number * nr_interleaves + interleaf + 1;
                            
                            phigold = 2.4;%137.51*2*pi/360;                              
                            phi = n*phigold;
                        
                            z = cos(sqrt(n/N) * pi/2);  
%                        %%     disp('phyllo classic ')
        if ssi==1
            if py_number == 0
                z = 1; 
                % disp('phyllo strictly SI first projection / interleaf')
            end
        end
%                        %% --------------------------------------------
                            
                            if mod(py_number,2)
                                z = z*-1;
                                phi = phi + pi;
                            end

%                             phi_all(py_number+1,interleaf+1) = rem(phi,pi*2);
%                             theta_all(py_number+1,interleaf+1) = (sqrt(n/N) * pi/2);
                phi_all(py_number+1,interleaf+1) = phi;
                theta_all(py_number+1,interleaf+1) = acos(z);
                            
                            Kpos(py_number+1, interleaf+1,1) = sin(phi)*sqrt(1-z^2);
                            Kpos(py_number+1, interleaf+1,2) = cos(phi)*sqrt(1-z^2);
                            Kpos(py_number+1, interleaf+1,3) = z;
%                             Kpos(:,py_number+1, interleaf+1,1) = sin(phi)*sqrt(1-z^2)*rad_pos;
%                             Kpos(:,py_number+1, interleaf+1,2) = cos(phi)*sqrt(1-z^2)*rad_pos;
%                             Kpos(:,py_number+1, interleaf+1,3) = z*rad_pos;
                            
                        end
end
                    
% phi_all = phi_all(:);
% theta_all = theta_all(:);