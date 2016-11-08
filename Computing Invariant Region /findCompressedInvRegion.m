function gCompressed = findCompressedInvRegion(F, g, A_c, B_c, Hu, delta_t)


    nFace = size(F,1) %number of faces of the region 

    nControlInputs = size(Hu,2)

    hu=ones(2*nControlInputs,1);


    uPoly = polytope(Hu,hu)
    Vu = extreme(uPoly)
    nControlExtremes = size(Vu, 1)



    gCompressed = g;
    FCompressed = F;

    %%
    while 1
        all_break_happened = 0;
        for face=1:nFace
            face
            compressionFactor = 0.99;
            compCoEff = 0.999;

            %for the current face find the maximum safety neighborhood region and 
            %construct the safe zone

            max_steps = 0;
            while 1
                max_steps = max_steps + 1;
                if max_steps > 200
                    display('cannot do it anymore!')
                    break;
                end

                orthFaceVector = F(face,:);

                gNeigh = [g; -compCoEff * g(face)];
                FNeigh = [F; -F(face,:)];
                neighPoly = polytope(FNeigh, gNeigh);
                V=extreme(neighPoly);

                nVertices = size(V,1);

                maxDer = -10000;

                for iV=1:nVertices
                    for iU=1:nControlExtremes

                        x_dot = A_c * V(iV,:)' + B_c * Vu(iU,:)';

                        derivitaveProjection = orthFaceVector * x_dot;
                        if derivitaveProjection > maxDer
                            maxDer = derivitaveProjection;
                            maxiV = iV;
                            maxiU = iU;
                        end
                    end
                end

    %             maxiV
    %             maxiU

                %now check if the maximum derivative is smaller than the time
                %required to pass through the neighborhood

                maxDistance = maxDer * delta_t;
                distance = abs(compCoEff * g(face) - g(face));

                if maxDistance <= distance
                    % if the neighborhood satisfies the conditions of derivatives
                    % then break and move to the next face of the neighborhood

                    gCompressed(face,:) = compCoEff * g(face,:);
                    all_break_happened = all_break_happened + 1;
                    break;
                end
    %             maxDistance
    %             distance
    %             maxDer

                compCoEff = compressionFactor * compCoEff;
            end

        end
        
        if all_break_happened == nFace
            break
    end
    
    %%
    
    

%     for face=1:nFace
%         face
%         compressionFactor = 0.99;
%         compCoEff = 0.999;
% 
%         %for the current face find the maximum safety neighborhood region and 
%         %construct the safe zone
%         
%         max_steps = 0;
%         while 1
%             max_steps = max_steps + 1;
%             if max_steps > 200
%                 display('cannot do it anymore!')
%                 break;
%             end
%             
%             orthFaceVector = F(face,:);
% 
%             gNeigh = [g; -compCoEff * g(face)];
%             FNeigh = [F; -F(face,:)];
%             neighPoly = polytope(FNeigh, gNeigh);
%             V=extreme(neighPoly);
% 
%             nVertices = size(V,1);
% 
%             maxDer = -10000;
% 
%             for iV=1:nVertices
%                 for iU=1:nControlExtremes
% 
%                     x_dot = A_c * V(iV,:)' + B_c * Vu(iU,:)';
% 
%                     derivitaveProjection = orthFaceVector * x_dot;
%                     if derivitaveProjection > maxDer
%                         maxDer = derivitaveProjection;
%                         maxiV = iV;
%                         maxiU = iU;
%                     end
%                 end
%             end
%             
% %             maxiV
% %             maxiU
% 
%             %now check if the maximum derivative is smaller than the time
%             %required to pass through the neighborhood
% 
%             maxDistance = maxDer * delta_t;
%             distance = abs(compCoEff * g(face) - g(face));
% 
%             if maxDistance <= distance
%                 % if the neighborhood satisfies the conditions of derivatives
%                 % then break and move to the next face of the neighborhood
% 
%                 gCompressed(face,:) = compCoEff * g(face,:);
% 
%                 break;
%             end
% %             maxDistance
% %             distance
% %             maxDer
% 
%             compCoEff = compressionFactor * compCoEff;
%         end
% 
%     end

end

