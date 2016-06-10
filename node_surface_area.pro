function mnf_make_data_array, nifti

    dim = nifti.dim
    type = nifti.datatype
    
;    dims = dim[1]
    
;    for i = 1, dim[0] do begin
;        dims = [ dims, dim[i+1] ne 0 ? dim[i+1] : 1 ]
;    endfor
    dims = dim[1:dim[0]]
    
    case type of
    
        0:  data_array = intarr(dims)       ;; unknown datatype
        1:  data_array = bytarr(dims)       ;; binary (1 bit/voxel)  (BT: unsure of IDLs format for this)
        2:  data_array = bytarr(dims)       ;; unsigned char (8 bits/voxel)
        4:  data_array = intarr(dims)       ;; signed short (16 bits/voxel)
        8:  data_array = lonarr(dims)       ;; signed int (32 bits/voxel)
        16: data_array = fltarr(dims)       ;; float (32 bits/voxel)
        32: data_array = complex_arr(dims)  ;; complex (64 bits/voxel)
        64: data_array = dblarr(dims)       ;; double (64 bits/voxel)
        256: data_array = bytarr(dims)      ;; signed char (8 bits)
        512: data_array = intarr(dims)      ;; unsigned short (16 bits)
        else: begin
            print, "Data type not supported."
            return, ptr_new()
        end
        
    endcase
    
    return, ptr_new(data_array, /no_copy)
;;    return, ptr_new(reform(data_array), /no_copy)
    
end

function mas_read_nifti, $
        nifti_filename  = nifti_filename, $
        default_dir = default_dir, $
        bval_filename  = bval_filename, $
        bvec_filename  = bvec_filename, $
        progress  = progress, $
        hdr_only  = hdr_only, $
        read_status=read_status
    
    catch, error_state
    if (error_state ne 0) then begin
        catch, /cancel
        void = dialog_message(['An error occurred while reading the NIFTI file:', $
                               !ERROR_STATE.msg], /center, /error)
        read_status = 0
        if (obj_valid(progressbar)) then progressbar->Destroy
        return, -1
    endif
    
    read_status = 0
    
    if (not keyword_set(nifti_filename)) then begin
    
        if (not keyword_set(default_dir)) then begin
            default_dir = project.current_path
        endif
        
        nifti_filename = dialog_pickfile(title='Select an .nii or .nii.gz file', $
            path=default_dir, filter='*.nii;*.nii.gz')
        if (nifti_filename eq '') then begin
            read_status = 0
            return, ptr_new()
        endif
        
    endif else begin
        if (not file_test(nifti_filename, /read)) then begin
            junk = dialog_message('Cannot read: '+nifti_filename, /error, /center)
            read_status = 0
            return, 0
        endif
    endelse
    
    if keyword_set(progress) and not keyword_set(hdr_only) then begin
        progressbar = obj_new('progressbar', title='MAS NIFTI READER', $
            text='Reading NIFTI file...', color='Red')
        progressbar->start
    endif
    
    if (not keyword_set(bval_filename)) then bval_filename = ''
    if (not keyword_set(bvec_filename)) then bvec_filename = ''
    
    openr, lun, nifti_filename, /get_lun, err=err, /compress
    
    if err ne 0 then begin
        read_status = 0
        return, 0
    endif
    
    read_status = 1
    
    nif_hdr = create_struct(name='MAS_NIFTI_HEADER')
    nifti   = create_struct(name='MAS_NIFTI_FILE')
    readu, lun, nif_hdr
    struct_assign, nif_hdr, nifti, /verbose
    nifti.nifti_filename = nifti_filename
    nifti.bval_filename = bval_filename
    nifti.bvec_filename = bvec_filename
    
    if (keyword_set(hdr_only)) then begin
        close, lun
        free_lun, lun
        return, nifti
    endif
    
    pdata = mnf_make_data_array(nifti)
    point_lun, lun, nifti.vox_offset

    blksize = long(nifti.dim[1])*long(nifti.dim[2])
    blk = (*pdata)[0:blksize-1] 
    nvoxels = n_elements(*pdata)
    for i = 0L, nvoxels-1, blksize do begin
        readu, lun, blk
        (*pdata)[i:i+blksize-1] = blk
        if (obj_valid(progressbar)) then begin
            progressbar->update, float(i)/float(nvoxels) * 100.0
        endif
    endfor
    ;;readu, lun, *pdata
    
    if (nifti.scl_slope ne 0) then begin
        *pdata = (*pdata * nifti.scl_slope) + nifti.scl_inter
    endif
    nifti.voxel_data = pdata
    
    close, lun
    free_lun, lun
    
    if (obj_valid(progressbar)) then begin
        progressbar->update, 100.0
        progressbar->destroy
    endif
    
    return, nifti
    
end

pro node_surface_area, node_file_name

    forward_function mas_read_nifti
    
    args = command_line_args()
    if (n_elements(args) eq 2) then begin
        node_file_name = args[0]
        save_file = args[1]
    endif else if (n_elements(node_file_name) eq 0) then begin
        node_file_name = dialog_pickfile(filter="*.nii;*.nii.gz", /read)
    endif
    
    if (node_file_name eq '') then return
    
    if (n_elements(save_file) eq 0) then begin
        save_file = dialog_pickfile(filter="*.csv", /write, default_extension='csv', /overwrite_prompt)
    endif
    
    nif = mas_read_nifti(nifti_filename=node_file_name, read_status=rs)
    
    if (rs eq 0) then begin
        junk = dialog_message("Unable to read nifti file.", /error, /center)
        return
    endif
    
    data = long(*nif.voxel_data)
    temp = lonarr(size(data, /dimensions))
    hist = histogram(data, locations=labels)
    
    voxel_tx = diag_matrix([nif.pixdim[1], nif.pixdim[2], nif.pixdim[3], 1.0])
    
    prog = widget_base(title="Progress", xoffset=300, yoffset=400, xpad=100, ypad=100)
    prog_b = widget_base(prog, xpad=50, ypad=50)
    txt = widget_label(prog_b, value="Please wait... calculating stats: 0000/"+string(n_elements(labels), format="(I04)"))
    widget_control, prog, /realize
    
    openw, lun, save_file, /get_lun
    printf, lun, "ROI Label,Voxel Count,Surface Area,Volume"
    for i = 1L, n_elements(labels)-1 do begin
        raw_ind = where(data eq labels[i], n_ind)
        if (n_ind eq 0) then begin
            ;;print, "Error: label"+string(labels[i], format='(I04)')+" has no voxels!"
            continue
        endif
        temp[raw_ind] = 100
        shade_volume, temp, 0, vertlist, polylist, /low;, /verbose
        
        vertlist = vert_t3d(vertlist, matrix=voxel_tx, /no_copy)
        
        ;o = obj_new('idlgrpolygon', vertlist, poly=polylist)
        ;xobjview, o
        
        printf, lun, strjoin([string(labels[i], format='(I04)'), $
                              n_elements(raw_ind), $
                              mesh_surfacearea(vertlist, polylist), $
                              n_elements(raw_ind)*product(nif.pixdim[1:3])], ",")
        temp[*] = 0
        
        if (i mod 20 eq 0) then begin
            widget_control, txt, set_value="Please wait... calculating stats: "+$
            string(i, format="(I0)")+"/"+string(n_elements(labels), format="(I0)")
        endif
        
    endfor
    widget_control, prog, /destroy
    close, lun
    free_lun, lun
    
    ptr_free, nif.voxel_data
            
end


pro mas_nifti_file__define

    struct = { MAS_NIFTI_FILE, inherits MAS_NIFTI_HEADER, $
        nifti_filename: '', $
        bval_filename: '', $
        bvec_filename: '', $
        voxel_data: ptr_new() }
        
end


pro mas_nifti_header__define

    struct = { MAS_NIFTI_HEADER, $
                sizeof_hdr: long(0), $
                data_type:  bytarr(10), $
                db_name: bytarr(18), $
                extents:  long(0), $
                session_error: 0, $
                regular: 0B, $
                dim_info: 0B, $
                dim: intarr(8), $
                intent_p1: 0.0, $
                intent_p2: 0.0, $
                intent_p3: 0.0, $
                intent_code: 0, $
                datatype: 0, $
                bitpix: 0, $
                slice_start: 0, $
                pixdim: fltarr(8), $
                vox_offset: 0.0, $
                scl_slope: 0.0, $
                scl_inter: 0.0, $
                slice_end: 0, $
                slice_code: 0B, $
                xyzt_units: 0B, $
                cal_max: 0.0, $
                cal_min: 0.0, $
                slice_duration: 0.0, $
                toffset: 0.0, $
                glmax: 0L, $
                glmin: 0L, $
                descrip: bytarr(80), $
                aux_file: bytarr(24), $
                qform_code: 0, $
                sform_code: 0, $
                quatern_b: 0.0, $
                quatern_c: 0.0, $
                quatern_d: 0.0, $
                qoffset_x: 0.0, $
                qoffset_y: 0.0, $
                qoffset_z: 0.0, $
                srow_x: fltarr(4), $
                srow_y: fltarr(4), $
                srow_z: fltarr(4), $
                intent_name: bytarr(16), $
                magic: bytarr(4) $
             }
        
end

    