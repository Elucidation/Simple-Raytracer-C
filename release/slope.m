function newarr = slope(arr)
newarr = zeros(size(arr));
for i = 1:length(arr)
    if (arr(i) == 0)
        if (i-1 > 0)
            left = arr(i-1);
            if left ~= 0
                newarr(i) = round(left / 2);
            end
        else
            left = 0;
        end
        if (i+1 > length(arr))
            right = 0;
        else
            right = arr(i+1);
            if right ~= 0
                newarr(i) = round(right / 2);
            end
        end
    else
        newarr(i) = arr(i);
    end
end
%return newarr;