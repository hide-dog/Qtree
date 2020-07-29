# ----------------------
# -- read             --
# ----------------------
function read_nodenum(skipnum)
    """ 
    xmax : 仮想セルも含めたnodeのxの数
    ymax : 仮想セルも含めたnodeのyの数
    """
    fff=[]
    open("grid/nodesnum", "r") do f
        fff=read(f,String)
    end 
    fff=split(fff,"\n",keepempty=false)
    num_nodes=length(fff)-skipnum

    for i in 1+skipnum:length(fff)
        fff[i]=replace(fff[i]," \r" => "")
    end

    temp = split(fff[2]," ")
    xmax = parse(Int64,temp[1]) 
    ymax = parse(Int64,temp[2]) 
    
    return xmax, ymax
end

function read_nodes(skipnum,xmax,ymax)
    """ 
    nodes[i][j][k]
    i : x点の番号
    j : y点の番号
    k=1 : 点のx座標
    k=2 : 点のy座標
    """

    fff=[]
    open("grid/nodes", "r") do f
        fff=read(f,String)
    end 
    fff=split(fff,"\n",keepempty=false)
    num_nodes=length(fff)-skipnum
    
    for i in 1+skipnum:length(fff)
        fff[i]=replace(fff[i]," \r" => "")
    end

    nodes = zeros(xmax,ymax,2)
    for i in 1:num_nodes
        temp=split(fff[i+skipnum]," ")

        xnum = parse(Int64,temp[1])
        ynum = parse(Int64,temp[2])
        nodes[xnum,ynum,1]=parse(Float64,temp[3])
        nodes[xnum,ynum,2]=parse(Float64,temp[4]) 
    end
    return nodes
end 

function read_nodes_vtk(skipnum)
    
    fff=[]
    open("grid/nodes_forvtk", "r") do f
        fff=read(f,String)
    end 
    fff=split(fff,"\n",keepempty=false)
    num_nodes=length(fff)-skipnum
    
    for i in 1+skipnum:length(fff)
        fff[i]=replace(fff[i]," \r" => "")
    end

    nodes=zeros(num_nodes,3)
    for i in 1:num_nodes
        temp=split(fff[i+skipnum]," ")

        # x = parse(Float64,temp[1])
        # y = parse(Float64,temp[2])
        # z = parse(Float64,temp[3])
        x = parse(Float64,temp[2])
        y = parse(Float64,temp[3])
        z = 0.0
        nodes[i,1] = x
        nodes[i,2] = y
        nodes[i,3] = z
    end
    return nodes
end 

function read_elements_vtk(skipnum)
    fff=[]
    open("grid/element_forvtk", "r") do f
        fff=read(f,String)
    end 
    fff=split(fff,"\n",keepempty=false)
    num_elements=length(fff)-skipnum
    
    for i in 1+skipnum:length(fff)
        fff[i]=replace(fff[i]," \r" => "")
    end

    elements = zeros(Int64,num_elements,4)
    for i in 1:num_elements
        temp=split(fff[i+skipnum]," ")
        
        elements[i,1] = parse(Int64,temp[2])
        elements[i,2] = parse(Int64,temp[3])
        elements[i,3] = parse(Int64,temp[4])
        elements[i,4] = parse(Int64,temp[5])
        
    end
    return elements
end 

function read_result(skipnum)
    fff=[]
    open("test333.dat", "r") do f
        fff=read(f,String)
    end 
    fff=split(fff,"\n",keepempty=false)
    num_point=length(fff)-skipnum
    
    for i in 1+skipnum:length(fff)
        fff[i]=replace(fff[i]," \r" => "")
    end

    readQ = zeros(num_point,5)
    for i in 1:num_point
        temp=split(fff[i+skipnum]," ")
        
        k = 1
        for j in 1:length(temp)
            if temp[j] != ""
                readQ[i,k] = parse(Float64,temp[j])
                k += 1
            end
        end

    end
    return readQ
end

function read_allgrid()
    skip=1
    xmax,ymax = read_nodenum(skip)
    nodes     = read_nodes(skip,xmax,ymax)
    nodes_vtk = read_nodes_vtk(skip)
    elements  = read_elements_vtk(skip)
    readQ     = read_result(skip)
    println("fin read grid")
    return xmax,ymax,nodes,nodes_vtk,elements,readQ
end
