function get_idlbl_path

; get the colormap directory from the idlbl_dir environmental variable
idlbl_dir = getenv('IDLBL_DIR')

; make sure the path ends with a `/'

if strlen(idlbl_dir) NE 0 then begin
    if strmid(idlbl_dir,strlen(idlbl_dir)-1,1) NE '/' then $
     idlbl_dir = idlbl_dir + '/'
endif

return, idlbl_dir
end
