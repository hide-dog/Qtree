using ProgressMeter

function main()
    PARAMDAT = "PARAMDAT.json"
    xmax,ymax,nodes,nodes_vtk,elements = read_allgrid()
    out_file,in_file,n_div,s_x,e_x,n_x,s_y,e_y,n_y = input_para(PARAMDAT)


    # 空間分割
    divnum = 2^n_div

    # モートン単位長
    unitx = abs(e_x - s_x) / divnum
    unity = abs(e_y - s_y) / divnum

    

    # 空間分割中心算出
    dx = e_x - s_x / n_x
    dy = e_y - s_y / n_y

    new_cell_center = zeros(n_x,n_y,2)
    new_morton_num = (n_x,n_y,2)

    for i in 1:n_x
        for j in 1:n_y
            # 左上の点のモートン空間番号
            ulx = Int((s_x + dx * (i-1) -s_x) / unitx) # 原点に合わせて番号を付けるため（-s_x）
            uly = Int((s_y + dy *  j    -s_y) / unity)

            # 右下の点のモートン空間番号
            lrx = Int((s_x + dx *  i    -s_x) / unitx) # 原点に合わせて番号を付けるため（-s_x）
            lry = Int((s_y + dy * (j-1) -s_y) / unity)

            # =-------------異常値をはじく事-----------

            morton_num[i,j,1],morton_num[i,j,2] = cal_morton(ulx,uly,lrx,lry)

            new_cell_center[i,j,1] = s_x +dx/2+dx*(i-1)
            new_cell_center[i,j,2] = s_y +dy/2+dy*(j-1)
        end
    end

    # セルを囲む四角形
    # これでモートン空間番号を算出
    over_rectangle = zeros(xmax-2,ymax-2,2,2)


    
    # セル中心算出
    cell_center = zeros(xmax-2,ymax-2,2)
    for i in 2:xmax-1
        for j in 2:ymax-1
            for k in 1:2
                cell_center[i,j,k] = 0.25*(nodes[i,j,k]+nodes[i,j+1,k]+nodes[i+1,j,k]+nodes[i+1,j+1,k])
            end
        end
    end

    

    
    

    




    
    
    prog = Progress(nt,1)
    for k in 1:nt
        next!(prog)

        
end

main()

