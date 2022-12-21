function index = bin_search(pos, array)
    low = 1;
    high = length(array);

    while (low < high)
        mid = floor(low + (high-low)/2);
        midVal = array(mid);

        if (midVal < pos)
            low = mid+1;
        else
            if (midVal >= pos)
                high = mid-1;
            end
        end
    end
    index= low; 
    
    if (array(index) > pos)
        index = index-1;
    end
end

