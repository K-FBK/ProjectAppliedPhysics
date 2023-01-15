classdef TTM
    %TTM- Transfer matrix method contains static methods
    %for constructing and transforming
    %wave-transfer and scettering matrices
    %
    % Karl-Fredrik Kylesten 2022

    methods(Static)
        function S = TransformSM(M)
            %MtoS Transforms wave transfer matrix M into Scattering matrix S
            S = (M(2,2))^(-1) * [M(1,1)*M(2,2)-M(2,1)*M(1,2), M(1,2); -M(2,1), 1];
        end

        function [t12,r12,t21,r21]=frensel_plane_TE(n1,n2,theta1,theta2)
            nominator_p = n1*cos(theta2)+n2*cos(theta1);
            t12 = 2*n1*cos(theta1)/nominator_p;
            t21 = 2*n2*cos(theta2)/nominator_p;
            r12 = (n2*cos(theta1)-n1*cos(theta2))/nominator_p;
            r21 = -r12;

        end

        function [t12,r12,t21,r21]=frensel_plane_TM(n1,n2,theta1,theta2)
           
            nominator_s = n1*cos(theta1)+n2*cos(theta2);
            t12 =  2*n1*cos(theta1)/nominator_s;
            t21 =  2*n2*cos(theta2)/nominator_s;
            r12 = (n1*cos(theta1)-n2*cos(theta2))/nominator_s;
            r21 = -r12;
            

        end

        function M = Mboundary(r12,t12,r21,t21)
            M = (t21^(-1))*[t12*t21-r12*r21, r21; -r12, 1];
        end

        function M = M_homogeneous(phi)
            %M_homogeneous creates a wave transfer matrix
            %for a homogeneous medium with phaseshift phi
            %phi = n*k_z*d*cos(theta)
            %d = thickness, k_0 = wavenumber, n = index
            
            M = [exp(1i*phi), 0; 0, exp(-1i*phi)]; 
        end
    end
end