# src path
src_path="C:\\Users\\Koike\\Desktop\\git\\Qtree\\src\\"

# main (変更しないこと)
src_read="read_grid.jl"
include(src_path*src_read)
src_read="read_para.jl"
include(src_path*src_read)
src_read="misc.jl"
include(src_path*src_read)
src_read="pre.jl"
include(src_path*src_read)


src_read="main_Qtree.jl"
include(src_path*src_read)


