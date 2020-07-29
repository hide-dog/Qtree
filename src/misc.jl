function ten_to_two(n)
    b = zeros(Int8,6)
    temp_b = zeros(Int8,6)

    ite=1
    while n > 0
        b[ite] = n % 2
        n = Int((n-b[ite])/2)
        ite = ite +1
    end
    
    #=
    for i in 1:6
        temp_b[i] = b[6+1-i]
    end
    =#
    
    return b
end

function two_to_ten(b)
    n = Int8(0)

    for i in 1:length(b)
        n = n + 2^(i-1)*b[i]
    end
    return n
end

function cal_morton(ulx,uly,lrx,lry) # 0<=s,tでreturn
    # 左上と右下の座標を整数で
    morton_belong = 0
    morton_num = 0

    upl_morton = zeros(Int8,6)
    lowr_morton = zeros(Int8,6)
    temp_morton = zeros(Int8,6)


    # 左上の点のモートン空間番号
    x_two = ten_to_two(ulx)
    y_two = ten_to_two(uly)
    
    for k in 1:3
        upl_morton[k*2-1] = x_two[k]
        upl_morton[k*2] = y_two[k]
    end

    # 右下の点のモートン空間番号
    x_two = ten_to_two(lrx)
    y_two = ten_to_two(lry)
    
    for k in 1:3
        lowr_morton[k*2-1] = x_two[k]
        lowr_morton[k*2] = y_two[k]
    end
    
    # 排他的論理和 xor
    for i in 1:6
        if upl_morton[i]+lowr_morton[i] == 2
            temp_morton[i] = 0
        else
            temp_morton[i] = upl_morton[i]+lowr_morton[i]
        end
    end
    
    # 所属空間のチェック
    if temp_morton[6] == 1 || temp_morton[5] == 1
        morton_belong = 0                             # ルート空間
    elseif temp_morton[4] == 1 || temp_morton[3] == 1
        morton_belong = 1                             # 親空間
    elseif temp_morton[2] == 1 || temp_morton[1] == 1
        morton_belong = 2                             # 子空間
    else
        morton_belong = 3                             # 孫空間
    end
    
    # 所属空間番号
    temp_morton = zeros(Int8,6)
    if morton_belong == 0                             # ルート空間
        morton_num = 0

    elseif morton_belong == 1                         # 親空間
        # 4右シフト
        temp_morton[1] = lowr_morton[5]
        temp_morton[2] = lowr_morton[6]
        
        morton_num = two_to_ten(temp_morton)

    elseif morton_belong == 2                         # 子空間
        # 2右シフト
        temp_morton[1] = lowr_morton[3]
        temp_morton[2] = lowr_morton[4]
        temp_morton[3] = lowr_morton[5]
        temp_morton[4] = lowr_morton[6]
        
        morton_num = two_to_ten(temp_morton)

    elseif morton_belong == 3                         # 孫空間
        temp_morton = copy(lowr_morton)
        morton_num = two_to_ten(temp_morton)
    end

    return morton_belong,morton_num
end

# ----------------------
# -- Inverse matrix   --
# ----------------------
function inverse_matrix(eu,ev)
# ----------------------
# -- A = ( eu[1]  ev[1])  --
# --     ( eu[2]  ev[2])  --
# ----------------------
    invA = zeros(2,2)
    detA = eu[1]*ev[2] - eu[2]*ev[1]

    invA[1,1] = (ev[2])/detA
    invA[1,2] = (-ev[1])/detA
    invA[2,1] = (-eu[2])/detA
    invA[2,2] = (eu[1])/detA
    return invA
end
    