using ProgressMeter

function main()
    PARAMDAT = "PARAMDAT.json"
    xmax,ymax,nodes,nodes_vtk,elements = read_allgrid()
    out_file,in_file,n_div,s_x,e_x,n_x,s_y,e_y,n_y = input_para(PARAMDAT)
    

    
    cell_center = zeros(xmax-2,ymax-2,2)
    for i in 2:xmax-1
        for j in 2:ymax-1
            for k in 1:2
                cell_center[i,j,k] = 0.25*(nodes[i,j,k]+nodes[i,j+1,k]+nodes[i+1,j,k]+nodes[i+1,j+1,k])
            end
        end
    end

    divnum = 4^n_div
    
    
    prog = Progress(nt,1)
    for k in 1:nt
        next!(prog)

        
end

main()