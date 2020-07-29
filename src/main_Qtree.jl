using ProgressMeter
using Printf

function main()
    PARAMDAT = "PARAMDAT.json"
    xmax,ymax,nodes,nodes_vtk,elements,readQ = read_allgrid()

    # println(readQ[1,:])
    


    out_file,in_file,n_div,s_x,e_x,n_x,s_y,e_y,n_y = input_para(PARAMDAT)

    lx = abs(e_x - s_x)
    ly = abs(e_y - s_y)

    main_pre(n_x,n_y,lx,ly)


    # 空間分割
    divnum = 2^n_div

    # モートン単位長
    unitx = abs(e_x - s_x) / divnum
    unity = abs(e_y - s_y) / divnum

    # 空間分割中心算出
    dx = (e_x - s_x) / n_x
    dy = (e_y - s_y) / n_y

    new_cell_center = zeros(n_x,n_y,2)
    new_morton_num = zeros(n_x,n_y,2)
    for i in 1:n_x
        for j in 1:n_y
            # 左上の点のモートン空間番号
            ulx = floor((s_x + dx * (i-1) -s_x) / unitx) # 原点に合わせて番号を付けるため（-s_x）
            uly = floor((s_y + dy *  j    -s_y) / unity)

            # 右下の点のモートン空間番号
            lrx = floor((s_x + dx *  i    -s_x) / unitx) # 原点に合わせて番号を付けるため（-s_x）
            lry = floor((s_y + dy * (j-1) -s_y) / unity)

            # =-------------異常値をはじく事-----------

            # 0<=s, 0<=t
            #new_morton_num[i,j,1],new_morton_num[i,j,2] = cal_morton(ulx,uly,lrx,lry)
            s,t = cal_morton(ulx,uly,lrx,lry)
            new_morton_num[i,j,1] = s
            new_morton_num[i,j,2] = t

            new_cell_center[i,j,1] = s_x +dx/2+dx*(i-1)
            new_cell_center[i,j,2] = s_y +dy/2+dy*(j-1)
        end
    end
    println(" fin kuukan bunkatu\n")

    # 登録管理構造体
    memori = 1000
    morton_tree = zeros(Int64,85,memori)
    morton_tree_num = ones(Int64,85)

    # セルを囲む四角形
    # これでモートン空間番号を算出
    cellnum = size(elements)[1]
    cell_center = zeros(cellnum,2)
    #morton_num = zeros(cellnum,2)
    
    for i in 1:cellnum
        x1 = nodes_vtk[elements[i,1],1]
        y1 = nodes_vtk[elements[i,1],2]
        x2 = nodes_vtk[elements[i,2],1]
        y2 = nodes_vtk[elements[i,2],2]
        x3 = nodes_vtk[elements[i,3],1]
        y3 = nodes_vtk[elements[i,3],2]
        x4 = nodes_vtk[elements[i,4],1]
        y4 = nodes_vtk[elements[i,4],2]
    
        temp = [x1,x2,x3,x4]
        maxx = maximum(temp)
        minx = minimum(temp)
        temp = [y1,y2,y3,y4]
        maxy = maximum(temp)
        miny = minimum(temp)

        # 左上の点のモートン空間番号
        ulx = floor((minx-s_x) / unitx) # 原点に合わせて番号を付けるため（-s_x）
        uly = floor((maxy-s_y) / unity)

        # 右下の点のモートン空間番号
        lrx = floor((maxx-s_x) / unitx) # 原点に合わせて番号を付けるため（-s_x）
        lry = floor((miny-s_y) / unity)

        # =-------------異常値をはじく事-----------
        
        #morton_num[i,1],morton_num[i,2] = cal_morton(ulx,uly,lrx,lry)
        
        # 0<=s, 0<=t
        s,t = cal_morton(ulx,uly,lrx,lry)

        num = Int64(t + (4^s-1)/3   +1)     # 線形空間番号

        s = morton_tree_num[num]
        morton_tree[num,s] = i
        morton_tree_num[num] += 1

        cell_center[i,1] = 0.25*(x1+x2+x3+x4)
        cell_center[i,2] = 0.25*(y1+y2+y3+y4)
    end
    println(" fin kuukan bunkatu2 \n")

    
    # 線形探索用
    # 両方の探査を行うため，煩雑
    search0 = zeros(Int64,85)
    search1 = zeros(Int64,4,22)
    search2 = zeros(Int64,16,7)
    search3 = zeros(Int64,64,4)
    for i in 1:85
        search0[i] = i
    end
    for i in 1:4
        search1[i,1] = 1　　　　　　　　　　　　　　　　# ルート空間
        search1[i,2] = i+1                           # 親空間
        for j in 1:4
            search1[i,j+2] = 6 + 4*(i-1) + (j-1)     # 子空間
        end
        for j in 1:16
            search1[i,j+6] = 22 + 16*(i-1) + (j-1)   # 孫空間
        end
    end
    for i in 1:16
        search2[i,1] = 1                             # ルート空間
        search2[i,2] = div(i-1,4) + 2                # 親空間
        search2[i,3] = i+5                           # 子空間
        for j in 1:4
            search2[i,j+3] = 22 + 4*(i-1) + (j-1)    # 孫空間
        end
    end
    for i in 1:64
        search3[i,1] = 1                             # ルート空間
        search3[i,2] = div(i-1,16) + 2               # 親空間
        search3[i,3] = div(i-1,4) + 6                # 子空間
        search3[i,4] = 22 + (i-1)                    # 孫空間
    end
    println(search1)
    println(search2)
    println(search3)

    # 分割したセル中心と近い点を三点計算
    close_point_num = zeros(n_x,n_y,3)
    distance = zeros(n_x,n_y,3)   
    distance[:,:,:] .= 100
    
    # 分割した空間がセルと当たるか否かの判定
    atari = zeros(Int64,n_x,n_y)
    
    prog = Progress(n_x,1)
    for i in 1:n_x
        next!(prog)
        for j in 1:n_y
            if new_morton_num[i,j,1] ==0      # 所属空間
                temp = copy(search0)
                loop = 85
            elseif new_morton_num[i,j,1] ==1      # 所属空間
                temp = copy(search1)
                loop = 88
            elseif new_morton_num[i,j,1] ==2      # 所属空間
                temp = copy(search2)
                loop = 112
            elseif new_morton_num[i,j,1] ==3      # 所属空間
                temp = copy(search3)
                loop = 256
            else 
                " stop program !!"
                throw(UndefVarError(:x))
            end

            # 空間事に探査
            for k in 1:loop
                for l in 1:memori
                    
                    if morton_tree[temp[k],l] == 0
                        break
                    end
                    bangou = morton_tree[temp[k],l]

                    x1 = new_cell_center[i,j,1]
                    y1 = new_cell_center[i,j,2]
                    x2 = cell_center[bangou,1]
                    y2 = cell_center[bangou,2]
                    temp_distance = ((x1-x2)^2 + (y1-y2)^2)^0.5

                    for ii in 1:3
                        if temp_distance < distance[i,j,ii]
                            for m in 1:(4-ii)-1
                                distance[i,j,4-m] = distance[i,j,3-m]
                                close_point_num[i,j,4-m] = close_point_num[i,j,3-m]
                            end
                            distance[i,j,ii] = temp_distance
                            close_point_num[i,j,ii] = morton_tree[temp[k],l]
                            break
                        end
                    end

                    # ---------- 過大評価 ---------------
                    cnum = morton_tree[temp[k],l]

                    x1 = nodes_vtk[elements[cnum,1],1]
                    y1 = nodes_vtk[elements[cnum,1],2]
                    x2 = nodes_vtk[elements[cnum,2],1]
                    y2 = nodes_vtk[elements[cnum,2],2]
                    x3 = nodes_vtk[elements[cnum,3],1]
                    y3 = nodes_vtk[elements[cnum,3],2]
                    x4 = nodes_vtk[elements[cnum,4],1]
                    y4 = nodes_vtk[elements[cnum,4],2]

                    temp_x = [x1,x2,x3,x4]
                    maxx = maximum(temp_x)
                    minx = minimum(temp_x)
                    temp_y = [y1,y2,y3,y4]
                    maxy = maximum(temp_y)
                    miny = minimum(temp_y)

                    x = new_cell_center[i,j,1]
                    y = new_cell_center[i,j,2]

                    if minx < x && x < maxx
                        if miny < y && y < maxy
                            atari[i,j] += 1
                        end
                    end
                end
            end
        end
    end

    println(close_point_num[23,25,:])


    println(" fin kuukan bunkatu3 \n")

    println("\n ---------------------------------- \n")
    println(" start interpolation ")
    println("\n ---------------------------------- \n")
    
    # 補間 loop
    Qcell = zeros(n_x,n_y,5)
    for i in 1:n_x
        for j in 1:n_y
            # cell中心ベクトル p
            p  = [new_cell_center[i,j,1], new_cell_center[i,j,2]]
            
            # 四点ベクトル p
            x0 = cell_center[Int64(close_point_num[i,j,1]),1]
            y0 = cell_center[Int64(close_point_num[i,j,1]),2]
            x1 = cell_center[Int64(close_point_num[i,j,2]),1]
            y1 = cell_center[Int64(close_point_num[i,j,2]),2]
            x2 = cell_center[Int64(close_point_num[i,j,3]),1]
            y2 = cell_center[Int64(close_point_num[i,j,3]),2]

            p0 = [x0, y0]
            p1 = [x1, y1]
            p2 = [x2, y2]
            
            # 局所座標ベクトル e
            eu = [p1[1] - p0[1], p1[2] - p0[2]]
            ev = [p2[1] - p0[1], p2[2] - p0[2]]
            inv_A = inverse_matrix(eu,ev)

            # 補間距離ベクトル u
            up = zeros(2)
            for l in 1:2
                up[l] = inv_A[l,1]*(p[1]-p0[1]) + inv_A[l,2]*(p[2]-p0[2])
            end
            for k in 1:5
                Qcell[i,j,k] = readQ[Int64(close_point_num[i,j,1]),k]
                            + up[1]*readQ[Int64(close_point_num[i,j,2]),k]
                            + up[2]*readQ[Int64(close_point_num[i,j,3]),k]
                #Qcell[i,j] = Qcell[i,j]*avoga
            end
        end
    end

    fff = "grid_morton/atari"
    open(fff,"w") do f
        write(f,"result:atari\n")
        for i in 1:n_x
            for j in 1:n_y
                a1 = string(atari[i,j])
                write(f, a1*"\n")
            end
        end
    end
    println("\nwrite "*fff)


    fff = "result_morton/"* out_file
    open(fff,"w") do f
        write(f,"result:rho[kg/m^3], u[m/s], v[m/s], p[Pa], T[K]\n")
        for i in 1:n_x
            for j in 1:n_y
    
                a1 = @sprintf("%8.8f", Qcell[i,j,1])
                a2 = @sprintf("%8.8f", Qcell[i,j,2])
                a3 = @sprintf("%8.8f", Qcell[i,j,3])
                a4 = @sprintf("%8.8f", Qcell[i,j,4])
                a5 = @sprintf("%8.8f", Qcell[i,j,5])

                #=
                a1 = string(Int64(close_point_num[i,j,1]))
                a2 = string(Int64(close_point_num[i,j,2]))
                a3 = string(Int64(close_point_num[i,j,3]))
                =#
    
                write(f, a1*" "*a2*" "*a3*" "*a4* " "*a5*"\n")
            end
        end
    end
    println("\nwrite "*fff)
        
end


# ----------------
main()
# ----------------
