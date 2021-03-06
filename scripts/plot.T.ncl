begin
;************************************************
; read in netCDF file
;************************************************

 nvar = "t2m"
 yr0 = 1979
 yr = 2003
 issn = 2
 version = (/"ssn","inter"/)
 do iv=0,1  ;dimsizes(version)-1

  a = addfile("../output/test/"+nvar+"."+version(iv)+".nc","r")
  ;************************************************
  ;if(any((/"qy","BetaM"/).eq.nvar(iv))) then
  lev = (/-10,10,20/)
  ;end if

  ;************************************************
    ssn=a->ssn
    season = ""+ssn(issn,0) +ssn(issn,1) +ssn(issn,2)

    varname=nvar+"_"+version(iv)
    var = a->$varname$
    printVarSummary(var)
    if (any(isnan_ieee(var))) then
     value = 1.e20
     replace_ieeenan (var, value, 0)
     var@_FillValue = value
    end if
  ;************************************************
  ; create default plot
  ;************************************************
    wks = gsn_open_wks("png","../output/test/"+nvar+"."+version(iv)+"."+season+yr)                  ; send graphics to PNG file

    res = True
    res@mpMinLatF            = -90      ; range to zoom in on
    res@mpMaxLatF            =  90.
    res@mpMinLonF            =  0.
    res@mpMaxLonF            = 360.

    res@cnFillOn             = True                  ; turn on color fill
    res@cnFillPalette        = "BlWhRe"              ; choose colormap


    res@cnLevelSelectionMode = "ManualLevels"	; manually set the contour levels with the following 3 resources
    res@cnMinLevelValF  = lev(0)			; set the minimum contour level
    res@cnMaxLevelValF  = lev(1)			; set the maximum contour level
    res@cnLevelSpacingF = lev(2)			; set the interval between contours

    res@lbLabelStride   = 2

    ;res@tiMainString         = var@long_name
    ;res@gsnCenterString = ""+ssn(issn,0) +ssn(issn,1) +ssn(issn,2)  + yr          ; plot center string
    res@tiMainString = season +yr

    print(nvar+version(iv))
    printVarSummary(var)

    print("plotting")
    plot = gsn_csm_contour_map(wks,var(yr-yr0,issn,:,:),res)         ; create a default plot
    print("end plotting")

end do
end
