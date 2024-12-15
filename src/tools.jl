
"""
    general_segm_approx(x, xs, ys, k)

    x: point to approximate
    xs: tuple/array of known X coordinates (in order, no duplecates)
    ys: tuple/array of known Y coordinates (in order, no duplecates)
    k_before: slope before the first point (to extrapolate)
    k_after: slope after the last point (to extrapolate)
"""
function general_segm_approx(x, xs, ys, k_before, k_after)
    # if need to exterpolate
    if x < xs[1]
        return k_before*(x-xs[1]) + ys[1]
    end
    if x > xs[end]
        return  k_after*(x-xs[end]) + ys[end]
    end

    # if need to interpolate
    for i in 1:length(xs)-1
        if xs[i] ≤ x ≤ xs[i+1]
            k_curr = (ys[i+1]-ys[i])/(xs[i+1]-xs[i])
            return k_curr*(x-xs[i]) + ys[i]
        end
    end
end

function general_segm_approx_after(x, xs, ys, k)
    k_before = (xs[2]-xs[1])/(ys[2]-ys[1])
    return general_segm_approx(x, xs, ys, k_before, k)
end

function general_segm_approx_symmetrical(x, xs, ys)
    k_before = (ys[2]-ys[1])/(xs[2]-xs[1])
    end_i = length(xs)
    k_after = (ys[end_i]-ys[end_i-1])/(xs[end_i]-xs[end_i-1])
    return general_segm_approx(x, xs, ys, k_before, k_after)
end

#########################################################################################

function sub_index(i_char)
    if i_char == '1' return "₁"
    elseif i_char == '2' return "₂"
    elseif i_char == '3' return "₃"
    elseif i_char == '4' return "₄"
    elseif i_char == '5' return "₅"
    elseif i_char == '6' return "₆"
    elseif i_char == '7' return "₇"
    elseif i_char == '8' return "₈"
    elseif i_char == '9' return "₉"
    elseif i_char == '0' return "₀"
    else return ""
    end
end

function sub_string(index::String)
    new_string = ""
    for i in eachindex(index)
        char = index[i]
        # println(char)
        # println(sub_index(char))
        new_string *= sub_index(char)
    end
    return new_string
end

function sub_string(index::Int)
    return sub_string(string(index))
end

# println("1234a")
# println(sub_string("1234a"))

#########################################################################################