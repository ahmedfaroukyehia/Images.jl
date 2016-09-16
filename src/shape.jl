# Hough Transform Implementation for detecting lines
function hough{T}(img::AbstractArray{T, 2}, precision::Number = 1.0, minPoints::Number = 500)
  thetaRange = -89:precision:90
  lenTheta = Float64(length(thetaRange))
  hyp = hypot(size(img.data,1), size(img.data,2))
  imgc = canny(img, 1, 0.995)

  houghMat = zeros(round(Integer,2*hyp)+1, round(Integer,lenTheta)+1)
  for j=1:size(imgc.data,2)
    for i=1:size(imgc.data,1)
      if imgc.data[i,j] > 0.9
        for th in thetaRange
          rho = i*cosd(th) + j*sind(th)
          houghMat[round(Integer,rho+hyp)+1, round(Integer,th/precision+(lenTheta/2))+1] += 1
        end
      end
    end
  end

  normHoughMat = houghMat/maximum(houghMat)
  params = Array{Float64}[]

  # coloring the lines in the original image
  for i=1:size(houghMat,1)
    for j=1:size(houghMat,2)
      if houghMat[i,j] > minPoints
        # calcualting the line parameters: y=ax+b
        a = tand(90+((j-(lenTheta/2)-1)*precision))
        b = (i-hyp-1)/cosd(90-((j-(lenTheta/2)-1)*precision))
        push!(params, [a b houghMat[i,j]])

        for j=1:size(img.data,2)
          for i=1:size(img.data,1)
            if abs(a*i + b -j) < 1
              img.data[i,j] = RGB(1,0,0)
            end
          end
        end
      end
    end
  end

  img, grayim(normHoughMat'), params
end
