begin
;************************************************
; read in netCDF file
;************************************************
 tnam = (/"t2m","ts"/)
 Tnam = (/"T2m","Ts"/)
 outvar = (/"corr","pval"/)
 invar = (/"r","pval"/)
 version = (/"","_dt"/)
 Version = (/"","Detrended "/)
 lat = (/ (/-90,-45/), (/45,90/) /)

  do it = 0,1
    do ih = 0,1
      do iv =0,1  ;version iv = 1 => detrended
         do ivar = 0,0 ;invar: corr or pval
           do is=0,3
            fin = "../output/corr"+version(iv)+"_"+tnam(it)+"_precip.nc"
            print(fin)
            a = addfile(fin,"r")

            ;************************************************
              ;if(any((/"qy","BetaM"/).eq.nvar(ivar))) then
              lev = (/-90, 90, 20/)
              ;end if

            ;************************************************
              ssn=a->ssn
              season = ""+ssn(is,0) +ssn(is,1) +ssn(is,2)

              varname=invar(ivar)
              var = a->$varname$
              printVarSummary(var)
              var=var*100.
              if (any(isnan_ieee(var))) then
               value = 1.e20
               replace_ieeenan (var, value, 0)
               var@_FillValue = value
              end if
            ;************************************************
            ; create default plot
            ;************************************************
              fout = "../output/"+outvar(ivar)+version(iv)+"_"+tnam(it)+"_precip."+season+".tropics"
              wks = gsn_open_wks("png",fout)                  ; send graphics to PNG file

              res = True

              res@mpMinLatF            = -30      ; range to zoom in on
              res@mpMaxLatF            =  30
              res@mpMinLonF            =  0.
              res@mpMaxLonF            = 360.

              res@mpShapeMode  = "FreeAspect"
              res@vpWidthF      = 0.8
              res@vpHeightF     = 0.2


              res@cnFillOn             = True                  ; turn on color fill
              res@cnFillPalette        = "BlAqGrWh2YeOrReVi22"              ; choose colormap


              res@cnLevelSelectionMode = "ManualLevels"	; manually set the contour levels with the following 3 resources
              res@cnMinLevelValF  = lev(0)			; set the minimum contour level
              res@cnMaxLevelValF  = lev(1)			; set the maximum contour level
              res@cnLevelSpacingF = lev(2)			; set the interval between contours

              res@lbLabelStride   = 2

              res@mpGeophysicalLineThicknessF = 2

              ;res@tiMainString         = var@long_name
              ;res@gsnCenterString = ""+ssn(issn,0) +ssn(issn,1) +ssn(issn,2)  + yr          ; plot center string
              res@tiMainString = Version(iv)+"ERAint "+ Tnam(it)+" vs GPCP precip" + ", "+season

              print("plotting")
              print(res)
              plot = gsn_csm_contour_map(wks,var(is,:,:),res)         ; create a default plot
              print("end plotting")
          end do
        end do
      end do
    end do
  end do
end
