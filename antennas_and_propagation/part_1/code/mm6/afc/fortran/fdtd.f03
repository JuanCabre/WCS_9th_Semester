! Fortran 2003 FDTD code
! v2.04 :: 07-05-2011
! (c) 2006- Ondrej Franek

! Variables shared with the interrupt subroutine
MODULE SHARED

IMPLICIT NONE

INTEGER, PARAMETER :: DOUBLE = KIND(0d0)

INTEGER(4) :: xmax, ymax, zmax, overall, freq_ff_size, freq_nf_size, resolution(2)
REAL(DOUBLE), DIMENSION(:), ALLOCATABLE :: out, energy
COMPLEX(DOUBLE), DIMENSION(:,:,:), ALLOCATABLE :: Eth, Eph
COMPLEX(DOUBLE), DIMENSION(:,:,:,:), ALLOCATABLE :: Ex_field, Ey_field, Ez_field
LOGICAL :: nfff_flag, field_flag, energy_flag

END MODULE SHARED



PROGRAM FDTD
        
USE SHARED

IMPLICIT NONE

INTEGER(1) :: byte

INTEGER(4) :: limits_1(3), limits_2(3), xmax_pml, ymax_pml, zmax_pml, pmax, pmax2, &
              nfff_flag_int, nfff_faces(6), field_flag_int, e_step, e_comp, length, Nt, &
              SourceSubs(3), SourceType, SourceDir, pad(3), padf(3), &
              i, j, k, l, id, jd, kd, t, p, m, n, x_wall, y_wall, z_wall, &
              x_nfff(2), y_nfff(2), z_nfff(2), &
              x_nfff_1(2), y_nfff_1(2), z_nfff_1(2), &
              nfff_size_x1, nfff_size_x2, &
              nfff_size_y1, nfff_size_y2, &
              nfff_size_z1, nfff_size_z2, &
              geom_center_2(3), delay, source_delay, source_count, probe_count, probe_number, &
              fieldsource_count, fieldsource_number, &
              drop_flag_int, drop2_flag_int, cache_flag_int, muaver_flag_int, &
              hardsource_flag_int, &
              time_start(8), time_end(8), partial_times(4), time_empty, &
              prec, overhead, pagesize, cachesize, assoc, cacheline, &
              geom_type, skip, round_point1(3), round_point2(3), round_point3(3), &
              materials, materials_prev, mat_number, lim1(3), lim2(3), component, &
              array_flags(4), array_size(3), max_vol, max_pl, &
              saved_element, subs_length, &
              media_i1, media_j1, media_k1, media_i1j1, media_j1k1, media_k1i1, media_here, &
              array_limits_1(3), array_limits_2(3), i2, j2, k2, direction, &
              tweaks_sizes(9), max_length, rise, fall
INTEGER(4), DIMENSION(:), ALLOCATABLE :: indx1, indy1, indz1, &
                                         indx2, indy2, indz2, &
                                         indx_pml,  indy_pml,  indz_pml,  &
                                         indx_pml1, indy_pml1, indz_pml1, &
                                         indx_pml2, indy_pml2, indz_pml2, &
                                         saved_column, comp, metal_lut, media_packed, &
                                         x_probe, y_probe, z_probe, probe_type, &
                                         fx1, fx2, fy1, fy2, fz1, fz2
INTEGER(4), DIMENSION(:,:), ALLOCATABLE :: saved_matrix, subs
INTEGER(4), DIMENSION(:,:,:), ALLOCATABLE :: media
                                         
INTEGER(8) :: filesize, filepos_geometry, filepos_object, filepos_start, filepos

REAL(DOUBLE) :: Dx, Dt, Dth, Dph, eta, c, pi, Rs, default_Rs, eps0, mu0, &
                phase_center(3), r_GP(3), &
                ph, sin_ph, cos_ph, th, sin_th, cos_th, &
                r_x, r_y, r_z, th_x, th_y, th_z, ph_x, ph_y, ph_z, &
                cos_temp1, sin_temp1, &
                cos_temp2, sin_temp2, &
                energy_now, energy_mult, energy_max, energy_percent, &
                Drop_pct, Drop2_pct, Drop_dB, Drop2_dB, &
                source_now, source_max, source_percent, &
                sig_max, kappa_max, a_max, m_s, m_k, m_a, &
                eps_r, sig, lut_actual(9), &
                point1(3), point2(3), point3(3), point4(3), p12(3), p32(3), &
                resistance, capacitance, normal(3), distance, tweak1, tweak2, &
                Cb0, Db0, s_temp(3), direction_double, fc, fs, mag, offset_dt
REAL(DOUBLE), DIMENSION(:), ALLOCATABLE :: Cax, Cay, Caz, &
                                           Cbx, Cby, Cbz, &
                                           Dbx, Dby, Dbz, &
                                           Cbx0, Cby0, Cbz0, &
                                           in, &
                                           b_e, b_h, c_e, c_h, &
                                           sig1, sig2, kappa1, kappa2, a1, a2, &
                                           indp1, indp2, &
                                           freq_ff, freq_nf, alpha, &
                                           x1_orig, x2_orig, &
                                           y1_orig, y2_orig, &
                                           z1_orig, z2_orig, &
                                           cos_x1, cos_y1, cos_z1, &
                                           cos_x2, cos_y2, cos_z2, &
                                           sin_x1, sin_y1, sin_z1, &
                                           sin_x2, sin_y2, sin_z2
REAL(DOUBLE), DIMENSION(:,:), ALLOCATABLE :: lut, lut_temp, probes
REAL(DOUBLE), DIMENSION(:,:,:), ALLOCATABLE :: Ex, Ey, Ez, &
                                               Hx, Hy, Hz, &
                                               Psi_Ex_y, Psi_Ex_z, &
                                               Psi_Ey_z, Psi_Ey_x, &
                                               Psi_Ez_x, Psi_Ez_y, &
                                               Psi_Hx_y, Psi_Hx_z, &
                                               Psi_Hy_z, Psi_Hy_x, &
                                               Psi_Hz_x, Psi_Hz_y

COMPLEX(DOUBLE) :: phase_factor, shift_GP, coef_GP, ex2, temp_H, temp_E
COMPLEX(DOUBLE), DIMENSION(:), ALLOCATABLE :: coef
COMPLEX(DOUBLE), DIMENSION(:,:), ALLOCATABLE :: Wff, Wnf
COMPLEX(DOUBLE), DIMENSION(:,:,:,:), ALLOCATABLE :: Ex_xy_dft, Ey_xy_dft, &
                                                    Ey_yz_dft, Ez_yz_dft, &
                                                    Ez_zx_dft, Ex_zx_dft, &
                                                    Hx_xy_dft, Hy_xy_dft, &
                                                    Hy_yz_dft, Hz_yz_dft, &
                                                    Hz_zx_dft, Hx_zx_dft, &
                                                    Hx_xy_dft2, Hy_xy_dft2, &
                                                    Hy_yz_dft2, Hz_yz_dft2, &
                                                    Hz_zx_dft2, Hx_zx_dft2
                                                    

LOGICAL :: pml_flag, drop_flag, drop2_flag, hardsource_flag, drop_cond, drop2_cond, &
           cache_flag, muaver_flag, probe_flag, lumpedsource_flag, fieldsource_flag
LOGICAL, DIMENSION(:,:), ALLOCATABLE :: mask_xy, mask_yz, mask_zx


TYPE FieldSourceType
	REAL(DOUBLE), DIMENSION(:,:,:), ALLOCATABLE :: magnitude, shift
	INTEGER(4) :: component, update
END TYPE FieldSourceType


TYPE(FieldSourceType), DIMENSION(:), ALLOCATABLE :: FieldSources


REAL :: time1, time2

INTEGER :: text_len, filename_len, cut_len, temp
CHARACTER :: eol(2) = [ACHAR(10), ACHAR(13)], chr, date, time, zone
CHARACTER, DIMENSION(:), ALLOCATABLE :: textfile, filename
CHARACTER(256) :: str_file

EXTERNAL INTERRUPT



DO

!! Reading a task from a file

OPEN (11, FILE='tasklist.txt', ACCESS='STREAM', POSITION='REWIND')

! First find the length of the file
i = 0
DO
    READ (11, END=100) chr
    i = i + 1
END DO
100 CONTINUE

text_len = i
IF (text_len==0) STOP 'No input file -- program terminated 1'

ALLOCATE( textfile(text_len) )

REWIND (11)
READ (11) textfile

DO i = 1,text_len
    IF (ANY(textfile(i)==eol)) EXIT
END DO
filename_len = i - 1

IF (filename_len==0) STOP 'No input file -- program terminated 2'

ALLOCATE( filename(filename_len) )

filename = textfile(1:filename_len)
textfile(1:filename_len) = eol(1)

DO i = 1,text_len
    IF (ALL(textfile(i)/=eol)) EXIT
END DO

REWIND (11)

IF (i <= text_len) THEN
    WRITE (11) textfile(i:)
ELSE
    WRITE (11) ''
END IF

CLOSE (11)

! Conversion to string
str_file = ''
DO i = 1,MIN(256,filename_len)
    str_file(i:i) = filename(i)
END DO

DEALLOCATE( textfile, filename )

PRINT *,     ''
PRINT '(A)', '**********************************************'  
PRINT *,     'Processing: ', TRIM(str_file)
PRINT '(A)', '**********************************************'  


CALL DATE_AND_TIME( date, time, zone, time_start )
CALL ELAPSED_TIME( time_empty )


! Reading file      
OPEN (11, FILE=TRIM(str_file), ACCESS='STREAM')
READ (11) limits_1, limits_2, pmax, nfff_flag_int, nfff_faces, field_flag_int, e_step, &
          drop_flag_int, drop2_flag_int, cache_flag_int, muaver_flag_int, &
          probe_count, fieldsource_count, Nt, x_nfff, y_nfff, z_nfff, &
          freq_ff_size, freq_nf_size, Dx, Dt, pi, default_Rs, eps0, mu0, fc, fs, mag, offset_dt

ALLOCATE( freq_ff(freq_ff_size), freq_nf(freq_nf_size), in(Nt), out(Nt) )
READ (11) in
          
! Various flags
pml_flag         = pmax              /= 0 ! PML present
nfff_flag        = nfff_flag_int     /= 0 ! Near to far field transform required
field_flag       = field_flag_int    /= 0 ! Near field DFT required
energy_flag      = e_step            /= 0 ! Energy required
drop_flag        = drop_flag_int     /= 0 ! Stop if total energy drops below threshold
drop2_flag       = drop2_flag_int    /= 0 ! Stop if source energy drops below threshold
cache_flag       = cache_flag_int    /= 0 ! Optimize for cache
muaver_flag      = muaver_flag_int   /= 0 ! Mu averaging
probe_flag       = probe_count       /= 0 ! Probes defined
fieldsource_flag = fieldsource_count /= 0 ! Field sources defined

IF (pml_flag)   READ (11) sig_max, kappa_max, a_max, m_s, m_k, m_a
IF (nfff_flag)  READ (11) resolution, eta, c, phase_center, freq_ff
IF (field_flag) READ (11) freq_nf
IF (drop_flag .OR. drop2_flag) READ (11) delay
IF (drop_flag)  READ (11) Drop_dB
IF (drop2_flag) READ (11) Drop2_dB, source_delay
IF (cache_flag) READ (11) prec, overhead, pagesize, cachesize, assoc, cacheline

xmax = ABS( limits_2(1) - limits_1(1) )
ymax = ABS( limits_2(2) - limits_1(2) )
zmax = ABS( limits_2(3) - limits_1(3) )

pmax2 = 2*pmax

xmax_pml = xmax + pmax2
ymax_pml = ymax + pmax2
zmax_pml = zmax + pmax2

max_length = xmax_pml*ymax_pml*zmax_pml*6/10
materials = 2 ! free space and PEC

ALLOCATE( media(xmax_pml,ymax_pml,zmax_pml) )
ALLOCATE( lut(max_length,9) )

media = 1
lut = 0
lut(1,:) = [1,1,1,1,1,1,0,0,0] ! free space
lut(2,:) = [0,0,0,1,1,1,0,0,0] ! PEC

INQUIRE (11, POS=filepos_geometry)
filepos_object = filepos_geometry



! First pass -- voxels
DO

    READ (11, POS=filepos_object) geom_type, skip
    filepos_object = filepos_object + skip + 8
    
    IF (skip == 0) EXIT
    
    SELECT CASE (geom_type)
    
        CASE (4) ! Brick
        
            READ (11) point1, point2, eps_r, sig
            
            round_point1 = NINT( point1, 4 )
            round_point2 = NINT( point2, 4 )
            
            mat_number = 0
            
            DO i = 1,materials
                IF (lut(i,1) == eps_r) THEN
                IF (lut(i,7) == sig)   THEN
                IF (lut(i,4) == 1)     THEN
                    mat_number = i
                    EXIT
                END IF
                END IF
                END IF
            END DO
        
            IF (mat_number == 0) THEN
                materials = materials + 1
                IF ( materials > max_length) STOP 'Error -- lookup table length exceeded'
                mat_number = materials
                lut(mat_number,1:3) = eps_r
                lut(mat_number,4:6) = 1
                lut(mat_number,7:9) = sig
            END IF

            lim1 = MAX( MIN( round_point1, round_point2 ), limits_1 ) - limits_1 + 1
            lim2 = MIN( MAX( round_point1, round_point2 ), limits_2 ) - limits_1
            
            IF ( ANY( media(lim1(1):lim2(1),lim1(2):lim2(2),lim1(3):lim2(3)) /= 1 ) ) &
                STOP 'Error -- objects overlapping'

            media(lim1(1):lim2(1),lim1(2):lim2(2),lim1(3):lim2(3)) = mat_number

                    
        CASE (8) ! External object
        
            READ (11) round_point1, array_flags
            IF ( array_flags(1) == 0 ) CYCLE ! no volume array
            
            READ (11) array_size
           
            lim1 = MAX( round_point1,              limits_1 ) - limits_1 + 1
            lim2 = MIN( round_point1 + array_size, limits_2 ) - limits_1

            round_point1 = round_point1 - limits_1 ! now it is in media coordinates
            
            max_vol = array_flags(1)
            IF ( materials + max_vol > max_length ) STOP 'Error -- lookup table length exceeded'

            INQUIRE (11, POS=filepos_start) ! get the starting position
            
            DO k = lim1(3),lim2(3)
            DO j = lim1(2),lim2(2)
            
                ! set the position in bytes
                filepos = filepos_start + (k-round_point1(3)-1)*array_size(1)*array_size(2) &
                                        + (j-round_point1(2)-1)*array_size(1) &
                                        + (lim1(1)-round_point1(1)-1)
                READ (11, POS=filepos)

                DO i = lim1(1),lim2(1)
                
                    READ (11) byte
                    IF ( byte /= 0 ) THEN
                        IF ( media(i,j,k) /= 1 ) &
                            STOP 'Error -- objects overlapping'
                        media(i,j,k) = INT( byte, 4 ) + materials
                    END IF
                
                END DO

            END DO
            END DO

            WHERE ( array_flags > 0 ) array_flags = 1
            
            ! find the first element after the arrays
            filepos = filepos_start + array_flags(1) * PRODUCT( array_size ) &
                                    + array_flags(2) * PRODUCT( array_size + [1,0,0] ) &
                                    + array_flags(3) * PRODUCT( array_size + [0,1,0] ) &
                                    + array_flags(4) * PRODUCT( array_size + [0,0,1] )
            
            READ (11, POS=filepos) lut( materials+1:materials+max_vol, [1,4,7] )
            
            lut(materials+1:materials+max_vol,[2,3,5,6,8,9]) = &
                lut(materials+1:materials+max_vol,[1,1,4,4,7,7])
            
            materials = materials + max_vol

    END SELECT    

END DO



! Detecting PEC bricks and setting them to material 2
DO i = 3,materials
    IF (lut(i,1) == 0) WHERE (media == i) media = 2
END DO



! Averaging
ALLOCATE( saved_matrix(xmax_pml,ymax_pml), saved_column(xmax_pml) )

saved_matrix = media(:,:,zmax_pml)
DO k = zmax_pml,1,-1

    saved_column = media(:,ymax_pml,k)
    DO j = ymax_pml,1,-1
    
        saved_element = media(xmax_pml,j,k)
        DO i = xmax_pml,1,-1
        
            media_here = media(i,j,k)

            IF (i == 1) THEN
                media_i1 = saved_element
                IF (j == 1) THEN
                    media_i1j1 = saved_column(xmax_pml)
                ELSE
                    media_i1j1 = media(xmax_pml,j-1,k)
                END IF
            ELSE
                media_i1 = media(i-1,j,k)
                IF (j == 1) THEN
                    media_i1j1 = saved_column(i-1)
                ELSE
                    media_i1j1 = media(i-1,j-1,k)
                END IF
            END IF
            
            IF (j == 1) THEN
                media_j1 = saved_column(i)
                IF (k == 1) THEN
                    media_j1k1 = saved_matrix(i,ymax_pml)
                ELSE
                    media_j1k1 = media(i,ymax_pml,k-1)
                END IF
            ELSE
                media_j1 = media(i,j-1,k)
                IF (k == 1) THEN
                    media_j1k1 = saved_matrix(i,j-1)
                ELSE
                    media_j1k1 = media(i,j-1,k-1)
                END IF
            END IF
            
            IF (k == 1) THEN
                media_k1 = saved_matrix(i,j)
                IF (i == 1) THEN
                    media_k1i1 = saved_matrix(xmax_pml,j)
                ELSE
                    media_k1i1 = saved_matrix(i-1,j)
                END IF
            ELSE
                media_k1 = media(i,j,k-1)
                IF (i == 1) THEN
                    media_k1i1 = media(xmax_pml,j,k-1)
                ELSE
                    media_k1i1 = media(i-1,j,k-1)
                END IF
            END IF

            ! eps_rx & sig_x
            IF ( media_here==2 .OR. media_j1==2 .OR. media_k1==2 .OR. media_j1k1==2 ) THEN
                lut_actual(1) = 0
                lut_actual(7) = 0
            ELSE
                lut_actual(1) = ( lut(media_here,1) + lut(media_j1,1) + lut(media_k1,1) + lut(media_j1k1,1) )/4
                lut_actual(7) = ( lut(media_here,7) + lut(media_j1,7) + lut(media_k1,7) + lut(media_j1k1,7) )/4
            END IF
                                
            ! eps_ry & sig_y
            IF ( media_here==2 .OR. media_k1==2 .OR. media_i1==2 .OR. media_k1i1==2 ) THEN
                lut_actual(2) = 0
                lut_actual(8) = 0
            ELSE
                lut_actual(2) = ( lut(media_here,2) + lut(media_k1,2) + lut(media_i1,2) + lut(media_k1i1,2) )/4
                lut_actual(8) = ( lut(media_here,8) + lut(media_k1,8) + lut(media_i1,8) + lut(media_k1i1,8) )/4
            END IF
                            
            ! eps_rz & sig_z
            IF ( media_here==2 .OR. media_i1==2 .OR. media_j1==2 .OR. media_i1j1==2 ) THEN
                lut_actual(3) = 0
                lut_actual(9) = 0
            ELSE
                lut_actual(3) = ( lut(media_here,3) + lut(media_i1,3) + lut(media_j1,3) + lut(media_i1j1,3) )/4
                lut_actual(9) = ( lut(media_here,9) + lut(media_i1,9) + lut(media_j1,9) + lut(media_i1j1,9) )/4
            END IF
                            
            ! mu_r x,y,z
            IF (muaver_flag) THEN
                lut_actual(4) = 2/( 1/lut(media_here,4) + 1/lut(media_i1,4) )
                lut_actual(5) = 2/( 1/lut(media_here,5) + 1/lut(media_j1,5) )
                lut_actual(6) = 2/( 1/lut(media_here,6) + 1/lut(media_k1,6) )
            ELSE
                lut_actual(4) = MAX( lut(media_here,4), lut(media_i1,4) )
                lut_actual(5) = MAX( lut(media_here,5), lut(media_j1,5) )
                lut_actual(6) = MAX( lut(media_here,6), lut(media_k1,6) )
            END IF    

            ! new type of material, preliminary...
            mat_number = materials + 1
            
            ! ...but try to find the same material among the previous lookup table entries
            DO
            
                ! is the neighbor in x direction the same?
                IF (i /= xmax_pml) THEN
                IF (ALL( lut_actual == lut(media(i+1,j,k),:) )) THEN
                    mat_number = media(i+1,j,k)
                    EXIT
                END IF
                END IF
                
                ! is the neighbor in y direction the same?
                IF (j /= ymax_pml) THEN
                IF (ALL( lut_actual == lut(media(i,j+1,k),:) )) THEN 
                    mat_number = media(i,j+1,k)
                    EXIT
                END IF
                END IF

                ! is the neighbor in z direction the same?
                IF (k /= zmax_pml) THEN
                IF (ALL( lut_actual == lut(media(i,j,k+1),:) )) THEN
                    mat_number = media(i,j,k+1)
                    EXIT
                END IF
                END IF
                
                ! ok, so let's check the whole table (worst case)
                DO n = 1,materials
                    IF ( lut_actual(1) == lut(n,1) .AND. &
                         lut_actual(2) == lut(n,2) .AND. &
                         lut_actual(3) == lut(n,3) .AND. &
                         lut_actual(4) == lut(n,4) .AND. &
                         lut_actual(5) == lut(n,5) .AND. &
                         lut_actual(6) == lut(n,6) .AND. &
                         lut_actual(7) == lut(n,7) .AND. &
                         lut_actual(8) == lut(n,8) .AND. &
                         lut_actual(9) == lut(n,9) ) THEN
                        mat_number = n
                        EXIT
                    END IF
                END DO
                
                EXIT
            
            END DO

            IF ( mat_number > materials ) THEN ! material is not, after all, present among the old ones
            
                materials = materials + 1
                IF ( materials > max_length ) STOP 'Error -- lookup table length exceeded'
                lut(mat_number,:) = lut_actual
                
            END IF

            media(i,j,k) = mat_number
            
        END DO
        
    END DO
    
END DO

DEALLOCATE( saved_matrix, saved_column )



! Return back to the beginning of the geometry
filepos_object = filepos_geometry

! Second pass -- resistor and capacitor
DO

    READ (11, POS=filepos_object) geom_type, skip
    filepos_object = filepos_object + skip + 8
    
    IF (skip == 0) EXIT
    
    SELECT CASE (geom_type)
    
        CASE (6) ! Resistor
        
            READ (11) point1, point2, resistance
            
            sig = ANINT( SUM( ABS( point1 - point2 ) ) ) / resistance / Dx
            
            CALL STAIRCASEWIRE( point1, point2 )

            materials_prev = materials ! previous number of materials
            
            DO i = 1,subs_length
            
                IF ( subs(i,1) > xmax_pml .OR. subs(i,2) > ymax_pml .OR. subs(i,3) > zmax_pml ) CYCLE

                lut_actual = lut(media(subs(i,1),subs(i,2),subs(i,3)),:)
                lut_actual(6 + ABS(comp(i))) = sig
                
                CALL SEARCH_LUT
                
                media(subs(i,1),subs(i,2),subs(i,3)) = mat_number
                
            END DO
            
            DEALLOCATE( subs, comp )
            

        CASE (7) ! Capacitor
        
            READ (11) point1, point2, capacitance
            
            eps_r = 1 + &
                    ( capacitance * 1e-9_DOUBLE ) * &
                    ANINT( SUM( ABS( point1 - point2 ) ) ) / eps0 / Dx
            
            CALL STAIRCASEWIRE( point1, point2 )

            materials_prev = materials ! previous number of materials
            
            DO i = 1,subs_length

                IF ( subs(i,1) > xmax_pml .OR. subs(i,2) > ymax_pml .OR. subs(i,3) > zmax_pml ) CYCLE

                lut_actual = lut(media(subs(i,1),subs(i,2),subs(i,3)),:)
                lut_actual(ABS(comp(i))) = eps_r
                
                CALL SEARCH_LUT

                media(subs(i,1),subs(i,2),subs(i,3)) = mat_number
                
            END DO
        
            DEALLOCATE( subs, comp )
            

    END SELECT    

END DO



SourceType = 0 ! 1 ... lumped, 2 ... field source, 0 ... source not specified

SourceDir = 0 ! reset of source direction

IF (probe_flag) THEN
    probe_number = 0 ! number of probe reset
    ALLOCATE( x_probe(probe_count), &
              y_probe(probe_count), &
              z_probe(probe_count), &
              probe_type(probe_count) )
END IF

IF (fieldsource_flag) THEN
    fieldsource_number = 0 ! number of probe reset
    ALLOCATE( FieldSources(fieldsource_count), &
    		  fx1(fieldsource_count), fx2(fieldsource_count), &
    		  fy1(fieldsource_count), fy2(fieldsource_count), &
    		  fz1(fieldsource_count), fz2(fieldsource_count) )
END IF

! Return back to the beginning of the geometry
filepos_object = filepos_geometry

! Third pass -- PEC objects and sources
DO 

    READ (11, POS=filepos_object) geom_type, skip
    filepos_object = filepos_object + skip + 8
    
    IF (skip == 0) EXIT
    
    SELECT CASE (geom_type)
    
        CASE (1) ! Wire
        
            READ (11) point1, point2
            
            CALL STAIRCASEWIRE( point1, point2 )

            materials_prev = materials ! previous number of materials
            
            DO i = 1,subs_length

                IF ( subs(i,1) > xmax_pml .OR. subs(i,2) > ymax_pml .OR. subs(i,3) > zmax_pml ) CYCLE

                lut_actual = lut(media(subs(i,1),subs(i,2),subs(i,3)),:)
                lut_actual(ABS(comp(i))) = 0
                
                CALL SEARCH_LUT
                
                media(subs(i,1),subs(i,2),subs(i,3)) = mat_number
                
            END DO
        
            DEALLOCATE( subs, comp )
            

        CASE (2,3) ! Triangle and Rectangle
        
            READ (11) point1, point2, point3
            
            SELECT CASE (geom_type)
            
                CASE (2) ! Triangle

                    array_limits_1 = NINT( MIN( point1, point2, point3 ), 4 )
                    array_limits_2 = NINT( MAX( point1, point2, point3 ), 4 )
                    
                CASE (3) ! Rectangle
            
                    IF (ALL(point1 == point2)) THEN ! Singular case
                        point3 = point1
                        point4 = point1
                    ELSE ! constructing the rectangle
                        p12 = point1 - point2
                        p32 = point3 - point2
                        point3 = point3 - p12 * ( DOT_PRODUCT( p12, p32 ) / DOT_PRODUCT( p12, p12 ) )
                        point4 = point3 + p12
                    END IF

                    array_limits_1 = NINT( MIN( point1, point2, point3, point4 ), 4 )
                    array_limits_2 = NINT( MAX( point1, point2, point3, point4 ), 4 )

            END SELECT
            
            array_size = array_limits_2 - array_limits_1 + 1
            
            ALLOCATE( mask_xy(array_size(1),array_size(2)), &
                      mask_yz(array_size(2),array_size(3)), &
                      mask_zx(array_size(3),array_size(1)) )
                      
            mask_xy = .FALSE.; mask_yz = .FALSE.; mask_zx = .FALSE.
            
            round_point1 = NINT( point1, 4 ) - array_limits_1 + 1
            
            CALL STAIRCASEWIRE( point1, point2 )
            DO i = 1,subs_length
                mask_xy(round_point1(1),round_point1(2)) = .TRUE.
                mask_yz(round_point1(2),round_point1(3)) = .TRUE.
                mask_zx(round_point1(3),round_point1(1)) = .TRUE.
                round_point1(ABS(comp(i))) = round_point1(ABS(comp(i))) + SIGN( 1, comp(i) )
            END DO
            DEALLOCATE( subs, comp )
            
            CALL STAIRCASEWIRE( point2, point3 )
            DO i = 1,subs_length
                mask_xy(round_point1(1),round_point1(2)) = .TRUE.
                mask_yz(round_point1(2),round_point1(3)) = .TRUE.
                mask_zx(round_point1(3),round_point1(1)) = .TRUE.
                round_point1(ABS(comp(i))) = round_point1(ABS(comp(i))) + SIGN( 1, comp(i) )
            END DO
            DEALLOCATE( subs, comp )
            
            SELECT CASE (geom_type)
            
                CASE (2) ! Triangle
                
                    CALL STAIRCASEWIRE( point3, point1 )
                    DO i = 1,subs_length
                        mask_xy(round_point1(1),round_point1(2)) = .TRUE.
                        mask_yz(round_point1(2),round_point1(3)) = .TRUE.
                        mask_zx(round_point1(3),round_point1(1)) = .TRUE.
                        round_point1(ABS(comp(i))) = round_point1(ABS(comp(i))) + SIGN( 1, comp(i) )
                    END DO
                    DEALLOCATE( subs, comp )
                    
                
                CASE (3) ! Rectangle
                
                    CALL STAIRCASEWIRE( point3, point4 )
                    DO i = 1,subs_length
                        mask_xy(round_point1(1),round_point1(2)) = .TRUE.
                        mask_yz(round_point1(2),round_point1(3)) = .TRUE.
                        mask_zx(round_point1(3),round_point1(1)) = .TRUE.
                        round_point1(ABS(comp(i))) = round_point1(ABS(comp(i))) + SIGN( 1, comp(i) )
                    END DO
                    DEALLOCATE( subs, comp )
                
                    CALL STAIRCASEWIRE( point4, point1 )
                    DO i = 1,subs_length
                        mask_xy(round_point1(1),round_point1(2)) = .TRUE.
                        mask_yz(round_point1(2),round_point1(3)) = .TRUE.
                        mask_zx(round_point1(3),round_point1(1)) = .TRUE.
                        round_point1(ABS(comp(i))) = round_point1(ABS(comp(i))) + SIGN( 1, comp(i) )
                    END DO
                    DEALLOCATE( subs, comp )

                                    
            END SELECT
            
            
            ! Fill the outlines
            
            DO j = 1,array_size(2)
                rise = 1
                fall = 0
                DO i = 1,array_size(1)   
                    IF ( mask_xy(i,j) ) THEN; rise = i; EXIT; END IF
                END DO
                DO i = array_size(1),1,-1
                    IF ( mask_xy(i,j) ) THEN; fall = i; EXIT; END IF
                END DO
                mask_xy(rise:fall,j) = .TRUE.                
            END DO
            
            DO j = 1,array_size(3)
                rise = 1
                fall = 0
                DO i = 1,array_size(2)   
                    IF ( mask_yz(i,j) ) THEN; rise = i; EXIT; END IF
                END DO
                DO i = array_size(2),1,-1
                    IF ( mask_yz(i,j) ) THEN; fall = i; EXIT; END IF
                END DO
                mask_yz(rise:fall,j) = .TRUE.                
            END DO
        
            DO j = 1,array_size(1)
                rise = 1
                fall = 0
                DO i = 1,array_size(3)   
                    IF ( mask_zx(i,j) ) THEN; rise = i; EXIT; END IF
                END DO
                DO i = array_size(3),1,-1
                    IF ( mask_zx(i,j) ) THEN; fall = i; EXIT; END IF
                END DO
                mask_zx(rise:fall,j) = .TRUE.                
            END DO

            
            ! Find patches
            
            DO j = 1,array_size(2)-1
            DO i = 1,array_size(1)-1
                ! If the four vertices of the patch are filled, mark
                mask_xy(i,j) = mask_xy(i  ,j  ) .AND. &
                               mask_xy(i+1,j  ) .AND. &
                               mask_xy(i  ,j+1) .AND. &
                               mask_xy(i+1,j+1)
            END DO
            END DO
        
            DO j = 1,array_size(3)-1
            DO i = 1,array_size(2)-1
                ! If the four vertices of the patch are filled, mark
                mask_yz(i,j) = mask_yz(i  ,j  ) .AND. &
                               mask_yz(i+1,j  ) .AND. &
                               mask_yz(i  ,j+1) .AND. &
                               mask_yz(i+1,j+1)
            END DO
            END DO
       
            DO j = 1,array_size(1)-1
            DO i = 1,array_size(3)-1
                ! If the four vertices of the patch are filled, mark
                mask_zx(i,j) = mask_zx(i  ,j  ) .AND. &
                               mask_zx(i+1,j  ) .AND. &
                               mask_zx(i  ,j+1) .AND. &
                               mask_zx(i+1,j+1)
            END DO
            END DO
            
            
            ! If the plane will not generate any patch, return
            IF ( .NOT. ANY( mask_xy ) .AND. .NOT. ANY( mask_yz ) .AND. .NOT. ANY( mask_zx ) ) CYCLE

            ! plane description:
            ! normal vector - decimated to single precision for compatibility
            ! between plates of different point order
            normal = CROSS( point1, point2 ) + CROSS( point2, point3 ) + CROSS( point3, point1 )
            normal = REAL( REAL( normal, 4 ), DOUBLE )
            
            IF ( normal(2) > 0 ) normal = -normal ! proper order of patches
            IF ( normal(3) > 0 ) normal = -normal ! for compatibility with wires

            ! distance from origin - shifted above zero
            distance = DOT_PRODUCT( normal, point1 - array_limits_1 + 0.5_DOUBLE )
            
            materials_prev = materials ! previous number of materials

            IF ( normal(3) /= 0 ) THEN
            
                DO j = 1,array_size(2)-1
                DO i = 1,array_size(1)-1
                
                    IF ( mask_xy(i,j) ) THEN    
                
                        k = CEILING( ( distance - normal(1)*i - normal(2)*j )/normal(3) , 4 )
                        IF ( normal(1)*i + normal(2)*j + normal(3)*k - distance > 0 ) k = k + 1
                        
                        ! media subscripts
                        i2 = i + array_limits_1(1) - limits_1(1)
                        j2 = j + array_limits_1(2) - limits_1(2)
                        k2 = k + array_limits_1(3) - limits_1(3)
                        
                        IF ( i2 > xmax_pml .OR. j2 > ymax_pml .OR. k2 > zmax_pml ) CYCLE

                        lut_actual = lut(media(i2,j2,k2),:)
                        IF ( lut_actual(1) /= 0 .OR. lut_actual(2) /= 0 ) THEN
                            lut_actual([1,2]) = 0
                            CALL SEARCH_LUT
                            media(i2,j2,k2) = mat_number
                        END IF

                        IF ( i2 < xmax_pml ) THEN
                            lut_actual = lut(media(i2+1,j2,k2),:)
                            IF ( lut_actual(2) /= 0 ) THEN
                                lut_actual(2) = 0
                                CALL SEARCH_LUT
                                media(i2+1,j2,k2) = mat_number
                            END IF
                        END IF
                        
                        IF ( j2 < ymax_pml ) THEN
                            lut_actual = lut(media(i2,j2+1,k2),:)
                            IF ( lut_actual(1) /= 0 ) THEN
                                lut_actual(1) = 0
                                CALL SEARCH_LUT
                                media(i2,j2+1,k2) = mat_number
                            END IF
                        END IF
                    
                    END IF

                END DO
                END DO
            
            END IF
            
            
            IF ( normal(1) /= 0 ) THEN
            
                DO k = 1,array_size(3)-1
                DO j = 1,array_size(2)-1
                
                    IF ( mask_yz(j,k) ) THEN    
                
                        IF ( normal(1) < 0 ) THEN
                            i = CEILING( ( distance - normal(2)*j - normal(3)*k )/normal(1) , 4 )
                            IF ( normal(1)*i + normal(2)*j + normal(3)*k - distance > 0 ) i = i + 1
                        ELSE ! ( normal(1) > 0 )
                            i = FLOOR( ( distance - normal(2)*j - normal(3)*k )/normal(1) , 4 ) + 1
                            IF ( normal(1)*i + normal(2)*j + normal(3)*k - distance < 0 ) i = i + 1
                        END IF
                        
                        ! media subscripts
                        i2 = i + array_limits_1(1) - limits_1(1)
                        j2 = j + array_limits_1(2) - limits_1(2)
                        k2 = k + array_limits_1(3) - limits_1(3)
                        
                        IF ( i2 > xmax_pml .OR. j2 > ymax_pml .OR. k2 > zmax_pml ) CYCLE

                        lut_actual = lut(media(i2,j2,k2),:)
                        IF ( lut_actual(2) /= 0 .OR. lut_actual(3) /= 0 ) THEN
                            lut_actual([2,3]) = 0
                            CALL SEARCH_LUT
                            media(i2,j2,k2) = mat_number
                        END IF

                        IF ( j2 < ymax_pml ) THEN
                            lut_actual = lut(media(i2,j2+1,k2),:)
                            IF ( lut_actual(3) /= 0 ) THEN
                                lut_actual(3) = 0
                                CALL SEARCH_LUT
                                media(i2,j2+1,k2) = mat_number
                            END IF
                        END IF
                        
                        IF ( k2 < zmax_pml ) THEN
                            lut_actual = lut(media(i2,j2,k2+1),:)
                            IF ( lut_actual(2) /= 0 ) THEN
                                lut_actual(2) = 0
                                CALL SEARCH_LUT
                                media(i2,j2,k2+1) = mat_number
                            END IF
                        END IF
                    
                    END IF

                END DO
                END DO
            
            END IF


            IF ( normal(2) /= 0 ) THEN
            
                DO k = 1,array_size(3)-1
                DO i = 1,array_size(1)-1
                
                    IF ( mask_zx(k,i) ) THEN
                
                        IF ( normal(2) < 0 ) THEN
                            j = CEILING( ( distance - normal(3)*k - normal(1)*i )/normal(2) , 4 )
                            IF ( normal(1)*i + normal(2)*j + normal(3)*k - distance > 0 ) j = j + 1
                        ELSE ! ( normal(2) > 0 )
                            j = FLOOR( ( distance - normal(3)*k - normal(1)*i )/normal(2) , 4 ) + 1
                            IF ( normal(1)*i + normal(2)*j + normal(3)*k - distance < 0 ) j = j + 1
                        END IF
                        
                        ! media subscripts
                        i2 = i + array_limits_1(1) - limits_1(1)
                        j2 = j + array_limits_1(2) - limits_1(2)
                        k2 = k + array_limits_1(3) - limits_1(3)
                        
                        IF ( i2 > xmax_pml .OR. j2 > ymax_pml .OR. k2 > zmax_pml ) CYCLE

                        lut_actual = lut(media(i2,j2,k2),:)
                        IF ( lut_actual(3) /= 0 .OR. lut_actual(1) /= 0 ) THEN
                            lut_actual([1,3]) = 0
                            CALL SEARCH_LUT
                            media(i2,j2,k2) = mat_number
                        END IF

                        IF ( k2 < zmax_pml ) THEN
                            lut_actual = lut(media(i2,j2,k2+1),:)
                            IF ( lut_actual(1) /= 0 ) THEN
                                lut_actual(1) = 0
                                CALL SEARCH_LUT
                                media(i2,j2,k2+1) = mat_number
                            END IF
                        END IF
                        
                        IF ( i2 < xmax_pml ) THEN
                            lut_actual = lut(media(i2+1,j2,k2),:)
                            IF ( lut_actual(3) /= 0 ) THEN
                                lut_actual(3) = 0
                                CALL SEARCH_LUT
                                media(i2+1,j2,k2) = mat_number
                            END IF
                        END IF
                    
                    END IF

                END DO
                END DO
            
            END IF

            
            DEALLOCATE( mask_xy, mask_yz, mask_zx )
            
            
        CASE (5) ! Source
        
            READ (11) point1, direction_double, resistance
            
            round_point1 = NINT( point1, 4 )
            
            SourceSubs = round_point1 - limits_1 + 1
            WHERE ( INT(direction_double) == [2,4,6] ) SourceSubs = SourceSubs - 1
            
            SourceDir = (INT(direction_double) + 1)/2
            IF ( SourceType == 0 ) SourceType = 1 ! lumped source only if not any field source
            
            Rs = resistance
            hardsource_flag = Rs == 0
            
            IF ( hardsource_flag ) THEN
                sig = 0
            ELSE
                sig = 1/(Rs*Dx)
            END IF
            
            lut_actual = lut(media(SourceSubs(1),SourceSubs(2),SourceSubs(3)),:)
            IF ( lut_actual(6 + SourceDir) /= sig ) THEN
            
                lut_actual(6 + SourceDir) = sig
                
                materials = materials + 1
                IF ( materials > max_length ) STOP 'Error -- lookup table length exceeded'
                
                lut(materials,:) = lut_actual
                media(SourceSubs(1),SourceSubs(2),SourceSubs(3)) = materials
                
            END IF
            
            
        CASE (9) ! Probe
        
            READ (11) round_point1, component
            
            round_point1 = round_point1 - limits_1 + 1
            probe_number = probe_number + 1
            
            x_probe(probe_number) = round_point1(1)
            y_probe(probe_number) = round_point1(2)
            z_probe(probe_number) = round_point1(3)
            probe_type(probe_number) = component
            
            
        CASE (8) ! External object
        
            READ (11) round_point1, array_flags
            
            max_vol = array_flags(1) ! length of the volume LUT
            max_pl  = MAXVAL( array_flags(2:4) ) ! length of the plate LUT
            
            WHERE ( array_flags > 0 ) array_flags = 1 ! normalization

            READ (11) array_size
           
            lim1 = MAX( round_point1,              limits_1 ) - limits_1 + 1
            lim2 = MIN( round_point1 + array_size, limits_2 ) - limits_1

            round_point1 = round_point1 - limits_1 ! now it is in media coordinates

            INQUIRE (11, POS=filepos_start) ! get the starting position

            ! Starting position for the plate arrays            
            filepos_start = filepos_start + array_flags(1) * PRODUCT( array_size )

            ! Jump the plate arrays and the volume LUTs
            filepos = filepos_start + array_flags(2) * PRODUCT( array_size + [1,0,0] ) &
                                    + array_flags(3) * PRODUCT( array_size + [0,1,0] ) &
                                    + array_flags(4) * PRODUCT( array_size + [0,0,1] ) &
                                    + array_flags(1) * ( max_vol * 3 * 8 )

            IF ( max_pl > 0 ) THEN 

                ! Read the metal LUT
                ALLOCATE( metal_lut(max_pl) )
                READ (11, POS=filepos) metal_lut
    
                materials_prev = materials ! previous number of materials
    
                
                ! Read the metal arrays
                IF ( array_flags(2) /= 0 ) THEN ! Plx array
    
                    DO k = lim1(3),lim2(3)
                    DO j = lim1(2),lim2(2)
        
                        ! set the position in bytes
                        filepos = filepos_start + (k-round_point1(3)-1)*(array_size(1)+1)*array_size(2) &
                                                + (j-round_point1(2)-1)*(array_size(1)+1) &
                                                + (lim1(1)-round_point1(1)-1)
                        READ (11, POS=filepos)
        
                        DO i = lim1(1),lim2(1)+1
                        
                            READ (11) byte
                            
                            IF ( byte == 0 ) CYCLE
                            IF ( metal_lut(byte) == 0 ) CYCLE
                            
                            lut_actual = lut(media(i,j,k),:)
                            
                            IF ( lut_actual(2) /= 0 .OR. lut_actual(3) /= 0 ) THEN
                                lut_actual([2,3]) = 0
                                CALL SEARCH_LUT
                                media(i,j,k) = mat_number
                            END IF
    
                            IF ( j < ymax_pml ) THEN
                                lut_actual = lut(media(i,j+1,k),:)
                                IF ( lut_actual(3) /= 0 ) THEN
                                    lut_actual(3) = 0
                                    CALL SEARCH_LUT
                                    media(i,j+1,k) = mat_number
                                END IF
                            END IF
                            
                            IF ( k < zmax_pml ) THEN
                                lut_actual = lut(media(i,j,k+1),:)
                                IF ( lut_actual(2) /= 0 ) THEN
                                    lut_actual(2) = 0
                                    CALL SEARCH_LUT
                                    media(i,j,k+1) = mat_number
                                END IF
                            END IF
                        
                        END DO
        
                    END DO
                    END DO
                    
                    filepos_start = filepos_start + PRODUCT( array_size + [1,0,0] )
                    
                END IF
                
                
                IF ( array_flags(3) /= 0 ) THEN ! Ply array
    
                    DO k = lim1(3),lim2(3)
                    DO j = lim1(2),lim2(2)+1
        
                        ! set the position in bytes
                        filepos = filepos_start + (k-round_point1(3)-1)*array_size(1)*(array_size(2)+1) &
                                                + (j-round_point1(2)-1)*array_size(1) &
                                                + (lim1(1)-round_point1(1)-1)
                        READ (11, POS=filepos)
        
                        DO i = lim1(1),lim2(1)
                        
                            READ (11) byte
                            
                            IF ( byte == 0 ) CYCLE
                            IF ( metal_lut(byte) == 0 ) CYCLE
                            
                            lut_actual = lut(media(i,j,k),:)
                            
                            IF ( lut_actual(3) /= 0 .OR. lut_actual(1) /= 0 ) THEN
                                lut_actual([1,3]) = 0
                                CALL SEARCH_LUT
                                media(i,j,k) = mat_number
                            END IF
    
                            IF ( k < zmax_pml ) THEN
                                lut_actual = lut(media(i,j,k+1),:)
                                IF ( lut_actual(1) /= 0 ) THEN
                                    lut_actual(1) = 0
                                    CALL SEARCH_LUT
                                    media(i,j,k+1) = mat_number
                                END IF
                            END IF
                            
                            IF ( i < xmax_pml ) THEN
                                lut_actual = lut(media(i+1,j,k),:)
                                IF ( lut_actual(3) /= 0 ) THEN
                                    lut_actual(3) = 0
                                    CALL SEARCH_LUT
                                    media(i+1,j,k) = mat_number
                                END IF
                            END IF
                        
                        END DO
        
                    END DO
                    END DO
                    
                    filepos_start = filepos_start + PRODUCT( array_size + [0,1,0] )
                    
                END IF
                
                
                IF ( array_flags(4) /= 0 ) THEN ! Plz array
    
                    DO k = lim1(3),lim2(3)+1
                    DO j = lim1(2),lim2(2)
        
                        ! set the position in bytes
                        filepos = filepos_start + (k-round_point1(3)-1)*array_size(1)*array_size(2) &
                                                + (j-round_point1(2)-1)*array_size(1) &
                                                + (lim1(1)-round_point1(1)-1)
                        READ (11, POS=filepos)
        
                        DO i = lim1(1),lim2(1)
                        
                            READ (11) byte
                            
                            IF ( byte == 0 ) CYCLE
                            IF ( metal_lut(byte) == 0 ) CYCLE
                            
                            lut_actual = lut(media(i,j,k),:)
                            
                            IF ( lut_actual(1) /= 0 .OR. lut_actual(2) /= 0 ) THEN
                                lut_actual([1,2]) = 0
                                CALL SEARCH_LUT
                                media(i,j,k) = mat_number
                            END IF
    
                            IF ( i < xmax_pml ) THEN
                                lut_actual = lut(media(i+1,j,k),:)
                                IF ( lut_actual(2) /= 0 ) THEN
                                    lut_actual(2) = 0
                                    CALL SEARCH_LUT
                                    media(i+1,j,k) = mat_number
                                END IF
                            END IF
                            
                            IF ( j < ymax_pml ) THEN
                                lut_actual = lut(media(i,j+1,k),:)
                                IF ( lut_actual(1) /= 0 ) THEN
                                    lut_actual(1) = 0
                                    CALL SEARCH_LUT
                                    media(i,j+1,k) = mat_number
                                END IF
                            END IF
                        
                        END DO
        
                    END DO
                    END DO
                    
                    filepos_start = filepos_start + PRODUCT( array_size + [0,0,1] )
                    
                END IF
    
    
                DEALLOCATE( metal_lut )
                
                filepos = filepos_start + array_flags(1) * ( max_vol * 3 * 8 ) &
                                        + max_pl * 4

                                
            END IF
            
            
            ! Source in external object
            READ (11, POS=filepos) direction
            
            IF ( direction > 0 ) THEN ! Source is present
            
                READ (11) round_point2
                
                SourceSubs = round_point2 + round_point1 + 1
                WHERE ( direction == [2,4,6] ) SourceSubs = SourceSubs - 1
                SourceDir = (direction + 1)/2
                IF ( SourceType == 0 ) SourceType = 1 ! lumped source only if not any field source
                
                Rs = default_Rs
                hardsource_flag = Rs == 0
                
                IF ( hardsource_flag ) THEN
                    sig = 0
                ELSE
                    sig = 1/(Rs*Dx)
                END IF
                
                lut_actual = lut(media(SourceSubs(1),SourceSubs(2),SourceSubs(3)),:)
                IF ( lut_actual(6 + SourceDir) /= sig ) THEN
                
                    lut_actual(6 + SourceDir) = sig
                    
                    materials = materials + 1
                    IF ( materials > max_length ) STOP 'Error -- lookup table length exceeded'
                    
                    lut(materials,:) = lut_actual
                    media(SourceSubs(1),SourceSubs(2),SourceSubs(3)) = materials
                    
                END IF
                
            END IF
        

        CASE (10) ! Field source
        
            READ (11) round_point1, round_point2
            
            lim1 = MIN( round_point1, round_point2 ) - limits_1 + 1
            lim2 = MAX( round_point1, round_point2 ) - limits_1 + 1

            ! check if out of bounds
            IF ( ANY( lim1 < 1 ) .OR. ANY( lim2 > [xmax,ymax,zmax] ) ) &
            	STOP 'Error -- field source out of bounds'
            	
            fieldsource_number = fieldsource_number + 1

            ! indices into the domain
            fx1(fieldsource_number) = lim1(1)
            fy1(fieldsource_number) = lim1(2)
            fz1(fieldsource_number) = lim1(3)
            fx2(fieldsource_number) = lim2(1)
            fy2(fieldsource_number) = lim2(2)
            fz2(fieldsource_number) = lim2(3)
            
            ! allocate
            ALLOCATE( FieldSources(fieldsource_number) % magnitude(lim2(1)-lim1(1)+1,lim2(2)-lim1(2)+1,lim2(3)-lim1(3)+1), &
            		  FieldSources(fieldsource_number) %     shift(lim2(1)-lim1(1)+1,lim2(2)-lim1(2)+1,lim2(3)-lim1(3)+1) )
            
            READ (11) FieldSources(fieldsource_number) % magnitude, &
                      FieldSources(fieldsource_number) % shift,     &
            		  FieldSources(fieldsource_number) % component, &
            		  FieldSources(fieldsource_number) % update
            		  
            SourceType = 2 ! Field source


    END SELECT    

END DO



IF ( SourceType == 0 ) STOP 'Error -- source not specified'

lumpedsource_flag = SourceType == 1



! Return back to the beginning of the geometry
filepos_object = filepos_geometry

! Fourth pass -- external object tweaks
DO 

    READ (11, POS=filepos_object) geom_type, skip
    filepos_object = filepos_object + skip + 8
    
    IF (skip == 0) EXIT
    
    SELECT CASE (geom_type)
    
        CASE (8) ! External object
        
            READ (11) round_point1, array_flags, array_size
            
            max_vol = array_flags(1) ! length of the volume LUT
            max_pl  = MAXVAL( array_flags(2:4) ) ! length of the plate LUT
            
            WHERE ( array_flags > 0 ) array_flags = 1 ! normalization

            round_point1 = round_point1 - limits_1 ! now it is in media coordinates

            INQUIRE (11, POS=filepos_start) ! get the starting position

            ! Jump to the source direction
            filepos = filepos_start + array_flags(1) * PRODUCT( array_size ) &
                                    + array_flags(2) * PRODUCT( array_size + [1,0,0] ) &
                                    + array_flags(3) * PRODUCT( array_size + [0,1,0] ) &
                                    + array_flags(4) * PRODUCT( array_size + [0,0,1] ) &
                                    + max_vol * 3 * 8 &
                                    + max_pl * 4

            ! Jump the source
            READ (11, POS=filepos) direction
            IF ( direction > 0 ) READ (11) round_point2
            
            ! Read the sizes of the tweaks
            READ (11) tweaks_sizes
            
            materials_prev = materials ! previous number of materials
            
            DO i = 1,9
                CALL ADD_TWEAKS(i)
            END DO

                    
    END SELECT
    
END DO



CLOSE (11)



! sigma at the source
IF ( lumpedsource_flag ) THEN
	mat_number = media(SourceSubs(1),SourceSubs(2),SourceSubs(3))
	sig = lut(mat_number,6 + SourceDir)
END IF

! Conversion of LUT to FDTD coefficients
Cb0 = (Dt / Dx) / eps0
Db0 = (Dt / Dx) / mu0
DO i = 1,materials

    WHERE ( lut(i,1:3) == 0 )
        lut_actual(1:3) = 1 ! Ca
        lut_actual(4:6) = 0 ! Cb
    ELSEWHERE
        s_temp = lut(i,7:9) * Dt / 2 / eps0 / lut(i,1:3)
        lut_actual(1:3) = ( 1 - s_temp ) / ( 1 + s_temp )   ! Ca
        lut_actual(4:6) = Cb0 / lut(i,1:3) / ( 1 + s_temp ) ! Cb
    END WHERE
    lut_actual(7:9) = Db0 / lut(i,4:6) ! Db
    lut(i,:) = lut_actual

END DO


IF ( lumpedsource_flag ) THEN

	! Input pulse calibration
	IF ( hardsource_flag ) THEN
	    in = -in / Dx
	ELSE
	    in = in * ( sig * lut(mat_number,3 + SourceDir) )
	END IF
	
END IF


! source energy termination condition cannot be used with field sources
IF ( fieldsource_flag ) drop2_flag = .FALSE.


IF (pml_flag) THEN

    ALLOCATE( b_e(pmax2), b_h(pmax2), c_e(pmax2), c_h(pmax2) )
    ALLOCATE(   sig1(pmax2),   sig2(pmax2), &
              kappa1(pmax2), kappa2(pmax2), &
                  a1(pmax2),     a2(pmax2), &
               indp1(pmax2),  indp2(pmax2)  )
              
    indp1 =   [ (i, i=0,pmax-1), (i, i=pmax  ,1,-1) ]                / REAL( pmax, DOUBLE )
    indp2 = ( [ (i, i=0,pmax-1), (i, i=pmax-1,0,-1) ] + 0.5_DOUBLE ) / REAL( pmax, DOUBLE )
    
    sig1 = sig_max * indp1 ** m_s
    sig2 = sig_max * indp2 ** m_s
    
    kappa1 = 1 + ( kappa_max - 1 ) * indp1 ** m_k
    kappa2 = 1 + ( kappa_max - 1 ) * indp2 ** m_k
    
    a1 = a_max * ( 1 - indp1 ) ** m_a
    a2 = a_max * ( 1 - indp2 ) ** m_a
    
    b_e = EXP( -( sig1 / kappa1 + a1 ) * Dt / eps0 )
    b_h = EXP( -( sig2 / kappa2 + a2 ) * Dt / eps0 )
    
    WHERE (sig1 == 0) ; c_e = 0
    ELSEWHERE         ; c_e = ( b_e - 1 ) / ( 1 + kappa1 * a1 * sig1**(-1) )
    END WHERE
                        c_h = ( b_h - 1 ) / ( 1 + kappa2 * a2 * sig2**(-1) )

    
    materials_prev = materials                  
                       
    ! x-oriented PML layer
    DO k = 1,zmax_pml
    DO j = 1,ymax_pml
    DO p = 1,pmax2
    
        i = xmax + p
        mat_number = media(i,j,k)
        lut_actual = lut(mat_number,:)
        
        lut_actual(4) = lut_actual(4) * kappa2(p)         ! Cbx
        lut_actual(5) = lut_actual(5) * ( 1 / kappa1(p) ) ! Cby
        lut_actual(6) = lut_actual(6) * ( 1 / kappa1(p) ) ! Cbz
        lut_actual(7) = lut_actual(7) * kappa1(p)         ! Dbx
        lut_actual(8) = lut_actual(8) * ( 1 / kappa2(p) ) ! Dby
        lut_actual(9) = lut_actual(9) * ( 1 / kappa2(p) ) ! Dbz

        ! new type of material, preliminary...
        mat_number = materials + 1
        
        ! ...but try to find the same material among the previous lookup table entries
        DO
        
            ! is the neighbor in y direction the same?
            IF (j > 1) THEN
            IF (ALL( lut_actual == lut(media(i,j-1,k),:) )) THEN 
                mat_number = media(i,j-1,k)
                EXIT
            END IF
            END IF

            ! is the neighbor in z direction the same?
            IF (k > 1) THEN
            IF (ALL( lut_actual == lut(media(i,j,k-1),:) )) THEN
                mat_number = media(i,j,k-1)
                EXIT
            END IF
            END IF
            
            ! ok, so let's check all the previous values in lut (worst case)
            DO n = materials_prev,materials
                IF ( lut_actual(1) == lut(n,1) .AND. &
                     lut_actual(2) == lut(n,2) .AND. &
                     lut_actual(3) == lut(n,3) .AND. &
                     lut_actual(4) == lut(n,4) .AND. &
                     lut_actual(5) == lut(n,5) .AND. &
                     lut_actual(6) == lut(n,6) .AND. &
                     lut_actual(7) == lut(n,7) .AND. &
                     lut_actual(8) == lut(n,8) .AND. &
                     lut_actual(9) == lut(n,9) ) THEN
                    mat_number = n
                    EXIT
                END IF
            END DO
            
            EXIT
        
        END DO

        IF ( mat_number > materials ) THEN ! material is not, after all, present among the old ones
        
            materials = materials + 1
            IF ( materials > max_length ) STOP 'Error -- lookup table length exceeded'
            lut(mat_number,:) = lut_actual
            
        END IF

        media(i,j,k) = mat_number

    END DO
    END DO
    END DO



    materials_prev = materials                  
    
    ! y-oriented PML layer
    DO k = 1,zmax_pml
    DO p = 1,pmax2
    DO i = 1,xmax_pml
    
        j = ymax + p
        mat_number = media(i,j,k)
        lut_actual = lut(mat_number,:)
        
        lut_actual(4) = lut_actual(4) * ( 1 / kappa1(p) ) ! Cbx
        lut_actual(5) = lut_actual(5) * kappa2(p)         ! Cby
        lut_actual(6) = lut_actual(6) * ( 1 / kappa1(p) ) ! Cbz
        lut_actual(7) = lut_actual(7) * ( 1 / kappa2(p) ) ! Dbx
        lut_actual(8) = lut_actual(8) * kappa1(p)         ! Dby
        lut_actual(9) = lut_actual(9) * ( 1 / kappa2(p) ) ! Dbz

        ! new type of material, preliminary...
        mat_number = materials + 1
        
        ! ...but try to find the same material among the previous lookup table entries
        DO
        
            ! is the neighbor in x direction the same?
            IF (i > 1) THEN
            IF (ALL( lut_actual == lut(media(i-1,j,k),:) )) THEN
                mat_number = media(i-1,j,k)
                EXIT
            END IF
            END IF

            ! is the neighbor in z direction the same?
            IF (k > 1) THEN
            IF (ALL( lut_actual == lut(media(i,j,k-1),:) )) THEN 
                mat_number = media(i,j,k-1)
                EXIT
            END IF
            END IF

            ! ok, so let's check all the previous values in lut (worst case)
            DO n = materials_prev,materials
                IF ( lut_actual(1) == lut(n,1) .AND. &
                     lut_actual(2) == lut(n,2) .AND. &
                     lut_actual(3) == lut(n,3) .AND. &
                     lut_actual(4) == lut(n,4) .AND. &
                     lut_actual(5) == lut(n,5) .AND. &
                     lut_actual(6) == lut(n,6) .AND. &
                     lut_actual(7) == lut(n,7) .AND. &
                     lut_actual(8) == lut(n,8) .AND. &
                     lut_actual(9) == lut(n,9) ) THEN
                    mat_number = n
                    EXIT
                END IF
            END DO
            
            EXIT
        
        END DO

        IF ( mat_number > materials ) THEN ! material is not, after all, present among the old ones
        
            materials = materials + 1
            IF ( materials > max_length ) STOP 'Error -- lookup table length exceeded'
            lut(mat_number,:) = lut_actual
            
        END IF

        media(i,j,k) = mat_number

    END DO
    END DO
    END DO



    materials_prev = materials                  

    ! z-oriented PML layer
    DO p = 1,pmax2
    DO j = 1,ymax_pml
    DO i = 1,xmax_pml
    
        k = zmax + p
        mat_number = media(i,j,k)
        lut_actual = lut(mat_number,:)
        
        lut_actual(4) = lut_actual(4) * ( 1 / kappa1(p) ) ! Cbx
        lut_actual(5) = lut_actual(5) * ( 1 / kappa1(p) ) ! Cby
        lut_actual(6) = lut_actual(6) * kappa2(p)         ! Cbz
        lut_actual(7) = lut_actual(7) * ( 1 / kappa2(p) ) ! Dbx
        lut_actual(8) = lut_actual(8) * ( 1 / kappa2(p) ) ! Dby
        lut_actual(9) = lut_actual(9) * kappa1(p)         ! Dbz

        ! new type of material, preliminary...
        mat_number = materials + 1
        
        ! ...but try to find the same material among the previous lookup table entries
        DO
        
            ! is the neighbor in x direction the same?
            IF (i > 1) THEN
            IF (ALL( lut_actual == lut(media(i-1,j,k),:) )) THEN 
                mat_number = media(i-1,j,k)
                EXIT
            END IF
            END IF

            ! is the neighbor in y direction the same?
            IF (j > 1) THEN
            IF (ALL( lut_actual == lut(media(i,j-1,k),:) )) THEN
                mat_number = media(i,j-1,k)
                EXIT
            END IF
            END IF
            
            ! ok, so let's check all the previous values in lut (worst case)
            DO n = materials_prev,materials
                IF ( lut_actual(1) == lut(n,1) .AND. &
                     lut_actual(2) == lut(n,2) .AND. &
                     lut_actual(3) == lut(n,3) .AND. &
                     lut_actual(4) == lut(n,4) .AND. &
                     lut_actual(5) == lut(n,5) .AND. &
                     lut_actual(6) == lut(n,6) .AND. &
                     lut_actual(7) == lut(n,7) .AND. &
                     lut_actual(8) == lut(n,8) .AND. &
                     lut_actual(9) == lut(n,9) ) THEN
                    mat_number = n
                    EXIT
                END IF
            END DO
            
            EXIT
        
        END DO

        IF ( mat_number > materials ) THEN ! material is not, after all, present among the old ones
        
            materials = materials + 1
            IF ( materials > max_length ) STOP 'Error -- lookup table length exceeded'
            lut(mat_number,:) = lut_actual
            
        END IF

        media(i,j,k) = mat_number

    END DO
    END DO
    END DO
    

    DEALLOCATE( sig1, sig2, kappa1, kappa2, a1, a2, indp1, indp2 )
    
END IF



ALLOCATE( lut_temp(materials,9) )
lut_temp = lut(1:materials,:)
DEALLOCATE( lut )

ALLOCATE( Cax(materials), Cay(materials), Caz(materials), &
          Cbx(materials), Cby(materials), Cbz(materials), &
          Dbx(materials), Dby(materials), Dbz(materials)  )

Cax = lut_temp(:,1)
Cay = lut_temp(:,2)
Caz = lut_temp(:,3)
Cbx = lut_temp(:,4)
Cby = lut_temp(:,5)
Cbz = lut_temp(:,6)
Dbx = lut_temp(:,7)
Dby = lut_temp(:,8)
Dbz = lut_temp(:,9)

DEALLOCATE( lut_temp )



pad  = [0,0,0]
padf = [0,0,0]
length = materials



! ALLOCATE( media_packed( SIZE(media) ) )
! media_packed = PACK( media, .TRUE. )
! OPEN (33, FILE=('output2.out'), STATUS='REPLACE', ACCESS='STREAM', FORM='UNFORMATTED', RECL=HUGE(0_8))
! WRITE (33) xmax_pml, ymax_pml, zmax_pml
! WRITE (33) Cax(media_packed), Cay(media_packed), Caz(media_packed), &
!            Cbx(media_packed), Cby(media_packed), Cbz(media_packed), &
!            Dbx(media_packed), Dby(media_packed), Dbz(media_packed)
! CLOSE (33)
! STOP 'output2.out written'


ALLOCATE( indx1(xmax_pml), indy1(ymax_pml), indz1(zmax_pml), &
          indx2(xmax_pml), indy2(ymax_pml), indz2(zmax_pml), &
          Ex(xmax_pml + pad(1),ymax_pml + pad(2),zmax_pml + pad(3)), &
          Ey(xmax_pml + pad(1),ymax_pml + pad(2),zmax_pml + pad(3)), &
          Ez(xmax_pml + pad(1),ymax_pml + pad(2),zmax_pml + pad(3)), &
          Hx(xmax_pml + pad(1),ymax_pml + pad(2),zmax_pml + pad(3)), &
          Hy(xmax_pml + pad(1),ymax_pml + pad(2),zmax_pml + pad(3)), &
          Hz(xmax_pml + pad(1),ymax_pml + pad(2),zmax_pml + pad(3))  )

IF (pml_flag) THEN
    ALLOCATE( indx_pml(pmax2),  indy_pml(pmax2),  indz_pml(pmax2),  &
              indx_pml1(pmax2), indy_pml1(pmax2), indz_pml1(pmax2), &
              indx_pml2(pmax2), indy_pml2(pmax2), indz_pml2(pmax2), &
              Psi_Ey_x(pmax2,ymax_pml,zmax_pml), &
              Psi_Ez_x(pmax2,ymax_pml,zmax_pml), &
              Psi_Ez_y(xmax_pml,pmax2,zmax_pml), &
              Psi_Ex_y(xmax_pml,pmax2,zmax_pml), &
              Psi_Ex_z(xmax_pml,ymax_pml,pmax2), &
              Psi_Ey_z(xmax_pml,ymax_pml,pmax2), &
              Psi_Hy_x(pmax2,ymax_pml,zmax_pml), &
              Psi_Hz_x(pmax2,ymax_pml,zmax_pml), &
              Psi_Hz_y(xmax_pml,pmax2,zmax_pml), &
              Psi_Hx_y(xmax_pml,pmax2,zmax_pml), &
              Psi_Hx_z(xmax_pml,ymax_pml,pmax2), &
              Psi_Hy_z(xmax_pml,ymax_pml,pmax2)  )
END IF

IF (nfff_flag) THEN

    x_nfff_1 = x_nfff + 1
    y_nfff_1 = y_nfff + 1
    z_nfff_1 = z_nfff + 1

    nfff_size_x2 = x_nfff(2) - x_nfff(1)
    nfff_size_y2 = y_nfff(2) - y_nfff(1)
    nfff_size_z2 = z_nfff(2) - z_nfff(1)

    nfff_size_x1 = nfff_size_x2 + 1
    nfff_size_y1 = nfff_size_y2 + 1
    nfff_size_z1 = nfff_size_z2 + 1

    ALLOCATE( Ex_xy_dft(  nfff_size_x2, nfff_size_y1, 2, freq_ff_size), &
              Ey_xy_dft(  nfff_size_x1, nfff_size_y2, 2, freq_ff_size), &
              Hx_xy_dft(  nfff_size_x1, nfff_size_y2, 2, freq_ff_size), &
              Hx_xy_dft2( nfff_size_x1, nfff_size_y2, 2, freq_ff_size), &
              Hy_xy_dft(  nfff_size_x2, nfff_size_y1, 2, freq_ff_size), &
              Hy_xy_dft2( nfff_size_x2, nfff_size_y1, 2, freq_ff_size), &
              Ey_yz_dft(  2, nfff_size_y2, nfff_size_z1, freq_ff_size), &
              Ez_yz_dft(  2, nfff_size_y1, nfff_size_z2, freq_ff_size), &
              Hy_yz_dft(  2, nfff_size_y1, nfff_size_z2, freq_ff_size), &
              Hy_yz_dft2( 2, nfff_size_y1, nfff_size_z2, freq_ff_size), &
              Hz_yz_dft(  2, nfff_size_y2, nfff_size_z1, freq_ff_size), &
              Hz_yz_dft2( 2, nfff_size_y2, nfff_size_z1, freq_ff_size), &
              Ez_zx_dft(  nfff_size_x1, 2, nfff_size_z2, freq_ff_size), &
              Ex_zx_dft(  nfff_size_x2, 2, nfff_size_z1, freq_ff_size), &
              Hz_zx_dft(  nfff_size_x2, 2, nfff_size_z1, freq_ff_size), &
              Hz_zx_dft2( nfff_size_x2, 2, nfff_size_z1, freq_ff_size), &
              Hx_zx_dft(  nfff_size_x1, 2, nfff_size_z2, freq_ff_size), &
              Hx_zx_dft2 (nfff_size_x1, 2, nfff_size_z2, freq_ff_size), &
              alpha(freq_ff_size), coef(freq_ff_size), &
              Wff(Nt,freq_ff_size) )
END IF

IF (field_flag) THEN
    ALLOCATE( Ex_field(xmax + padf(1),ymax + padf(2),zmax + padf(3),freq_nf_size), &
              Ey_field(xmax + padf(1),ymax + padf(2),zmax + padf(3),freq_nf_size), &
              Ez_field(xmax + padf(1),ymax + padf(2),zmax + padf(3),freq_nf_size), &
              Wnf(Nt,freq_nf_size) )
END IF

IF (energy_flag) THEN
    ALLOCATE( energy(Nt), Cbx0(length), Cby0(length), Cbz0(length) )
END IF

IF (probe_flag) THEN
    ALLOCATE( probes(Nt,probe_count) )
END IF

! auxiliary indices
indx1 = [ xmax_pml, (i, i=1,xmax_pml-1) ]
indy1 = [ ymax_pml, (i, i=1,ymax_pml-1) ]
indz1 = [ zmax_pml, (i, i=1,zmax_pml-1) ]

indx2 = [ (i, i=2,xmax_pml), 1 ]
indy2 = [ (i, i=2,ymax_pml), 1 ]
indz2 = [ (i, i=2,zmax_pml), 1 ]

IF (pml_flag) THEN
    indx_pml = [ (i, i=xmax+1,xmax_pml) ]
    indy_pml = [ (i, i=ymax+1,ymax_pml) ]
    indz_pml = [ (i, i=zmax+1,zmax_pml) ]

    indx_pml1 = [ (i, i=xmax,xmax_pml-1) ]
    indy_pml1 = [ (i, i=ymax,ymax_pml-1) ]
    indz_pml1 = [ (i, i=zmax,zmax_pml-1) ]

    indx_pml2 = [ (i, i=xmax+2,xmax_pml), 1 ]
    indy_pml2 = [ (i, i=ymax+2,ymax_pml), 1 ]
    indz_pml2 = [ (i, i=zmax+2,zmax_pml), 1 ]

    x_wall = xmax + 1 + pmax
    y_wall = ymax + 1 + pmax
    z_wall = zmax + 1 + pmax
ELSE
    x_wall = 1
    y_wall = 1
    z_wall = 1
END IF


! Array initialization
Ex = 0; Ey = 0; Ez = 0; Hx = 0; Hy = 0; Hz = 0; out = 0

IF (pml_flag) THEN
    Psi_Ey_x = 0; Psi_Ez_x = 0; Psi_Ez_y = 0
    Psi_Ex_y = 0; Psi_Ex_z = 0; Psi_Ey_z = 0
    Psi_Hy_x = 0; Psi_Hz_x = 0; Psi_Hz_y = 0
    Psi_Hx_y = 0; Psi_Hx_z = 0; Psi_Hy_z = 0
END IF

IF (nfff_flag) THEN 
    Ex_xy_dft = 0; Ey_xy_dft = 0
    Ey_yz_dft = 0; Ez_yz_dft = 0
    Ez_zx_dft = 0; Ex_zx_dft = 0 
    Hx_xy_dft = 0; Hy_xy_dft = 0
    Hy_yz_dft = 0; Hz_yz_dft = 0
    Hz_zx_dft = 0; Hx_zx_dft = 0
    Hx_xy_dft2 = 0; Hy_xy_dft2 = 0
    Hy_yz_dft2 = 0; Hz_yz_dft2 = 0
    Hz_zx_dft2 = 0; Hx_zx_dft2 = 0
END IF

IF (field_flag) THEN
    Ex_field = 0; Ey_field = 0; Ez_field = 0
END IF

IF (energy_flag) THEN
    energy = 0
END IF

IF (probe_flag) THEN
    probes = 0
END IF

IF (nfff_flag) THEN
    DO p = 1,freq_ff_size
    DO t = 1,Nt
        Wff(t,p) = EXP( CMPLX(0.0, -2.0*pi*((t - 1)*freq_ff(p))*Dt, DOUBLE) )
    END DO
    END DO
END IF

IF (field_flag) THEN
    DO p = 1,freq_nf_size
    DO t = 1,Nt
        Wnf(t,p) = EXP( CMPLX(0.0, -2.0*pi*((t - 1)*freq_nf(p))*Dt, DOUBLE) )
    END DO
    END DO
END IF

IF (energy_flag) THEN
    Cbx0 = Cbx*(2/(1 + Cax)); WHERE (Cbx0 == 0) Cbx0 = 1
    Cby0 = Cby*(2/(1 + Cay)); WHERE (Cby0 == 0) Cby0 = 1
    Cbz0 = Cbz*(2/(1 + Caz)); WHERE (Cbz0 == 0) Cbz0 = 1
    energy_mult = Dt*(Dx**2)/2
    e_comp = e_step;
    energy_max = 0;
END IF

IF (drop_flag) THEN
    Drop_pct = 100 * 10**(Drop_dB/10)
END IF

IF (drop2_flag) THEN
    Drop2_pct = 100 * 10**(Drop2_dB/10)
    source_now = 0
    source_max = 0
    source_count = 0
    source_percent = 0
END IF

energy_now = 0
energy_percent = 0

drop_cond = .TRUE.
drop2_cond = .TRUE.



CALL SIGNAL(2,INTERRUPT,temp)

CALL CPU_TIME( time1 )
CALL ELAPSED_TIME( partial_times(1) )



! FDTD loop
DO t = 1,Nt

!!! E component update !!!

! regular region
DO k=1,zmax_pml
DO j=1,ymax_pml
DO i=1,xmax_pml
    Ex(i,j,k) = Cax(media(i,j,k))*Ex(i,j,k) + &
                Cbx(media(i,j,k))*( Hz(i,j,k) - Hz(i,indy1(j),k) &
                                  - Hy(i,j,k) + Hy(i,j,indz1(k)) )
    Ey(i,j,k) = Cay(media(i,j,k))*Ey(i,j,k) + &
                Cby(media(i,j,k))*( Hx(i,j,k) - Hx(i,j,indz1(k)) &
                                  - Hz(i,j,k) + Hz(indx1(i),j,k) )
    Ez(i,j,k) = Caz(media(i,j,k))*Ez(i,j,k) + &
                Cbz(media(i,j,k))*( Hy(i,j,k) - Hy(indx1(i),j,k) &
                                  - Hx(i,j,k) + Hx(i,indy1(j),k) )
END DO
END DO
END DO

! PML region
IF (pml_flag) THEN

DO k=1,zmax_pml
DO j=1,ymax_pml
DO i=1,pmax2
    Psi_Ey_x(i,j,k) = b_e(i)*Psi_Ey_x(i,j,k) &
        + c_e(i)*( Hz(indx_pml(i),j,k) - Hz(indx_pml1(i),j,k) )
    Ey(indx_pml(i),j,k) = Ey(indx_pml(i),j,k) &
        - Cby(media(indx_pml(i),j,k))*Psi_Ey_x(i,j,k)
    Psi_Ez_x(i,j,k) = b_e(i)*Psi_Ez_x(i,j,k) &
        + c_e(i)*( Hy(indx_pml(i),j,k) - Hy(indx_pml1(i),j,k) )
    Ez(indx_pml(i),j,k) = Ez(indx_pml(i),j,k) &
        + Cbz(media(indx_pml(i),j,k))*Psi_Ez_x(i,j,k)
END DO
END DO
END DO

DO k=1,zmax_pml
DO j=1,pmax2
DO i=1,xmax_pml
    Psi_Ez_y(i,j,k) = b_e(j)*Psi_Ez_y(i,j,k) &
        + c_e(j)*( Hx(i,indy_pml(j),k) - Hx(i,indy_pml1(j),k) )
    Ez(i,indy_pml(j),k) = Ez(i,indy_pml(j),k) &
        - Cbz(media(i,indy_pml(j),k))*Psi_Ez_y(i,j,k)
    Psi_Ex_y(i,j,k) = b_e(j)*Psi_Ex_y(i,j,k) &
        + c_e(j)*( Hz(i,indy_pml(j),k) - Hz(i,indy_pml1(j),k) )
    Ex(i,indy_pml(j),k) = Ex(i,indy_pml(j),k) &
        + Cbx(media(i,indy_pml(j),k))*Psi_Ex_y(i,j,k)
END DO
END DO
END DO

DO k=1,pmax2
DO j=1,ymax_pml
DO i=1,xmax_pml
    Psi_Ex_z(i,j,k) = b_e(k)*Psi_Ex_z(i,j,k) &
        + c_e(k)*( Hy(i,j,indz_pml(k)) - Hy(i,j,indz_pml1(k)) )
    Ex(i,j,indz_pml(k)) = Ex(i,j,indz_pml(k)) &
        - Cbx(media(i,j,indz_pml(k)))*Psi_Ex_z(i,j,k)
    Psi_Ey_z(i,j,k) = b_e(k)*Psi_Ey_z(i,j,k) &
        + c_e(k)*( Hx(i,j,indz_pml(k)) - Hx(i,j,indz_pml1(k)) )
    Ey(i,j,indz_pml(k)) = Ey(i,j,indz_pml(k)) &
        + Cby(media(i,j,indz_pml(k)))*Psi_Ey_z(i,j,k)
END DO
END DO
END DO

END IF

! PEC walls
Ey(x_wall,:,:) = 0
Ez(x_wall,:,:) = 0
Ez(:,y_wall,:) = 0
Ex(:,y_wall,:) = 0
Ex(:,:,z_wall) = 0
Ey(:,:,z_wall) = 0



! Field sources
IF (fieldsource_flag) THEN

	DO l = 1,fieldsource_count
	
		SELECT CASE ( FieldSources(l) % component )
		
			CASE (1) ! Ex
				IF ( FieldSources(l) % update == 1 ) Ex(fx1(l):fx2(l),fy1(l):fy2(l),fz1(l):fz2(l)) = 0
				Ex(fx1(l):fx2(l),fy1(l):fy2(l),fz1(l):fz2(l)) = &
				Ex(fx1(l):fx2(l),fy1(l):fy2(l),fz1(l):fz2(l)) + &
				FieldSources(l) % magnitude * GAUSS( t + FieldSources(l) % shift )
			
			CASE (2) ! Ey
				IF ( FieldSources(l) % update == 1 ) Ey(fx1(l):fx2(l),fy1(l):fy2(l),fz1(l):fz2(l)) = 0
				Ey(fx1(l):fx2(l),fy1(l):fy2(l),fz1(l):fz2(l)) = &
				Ey(fx1(l):fx2(l),fy1(l):fy2(l),fz1(l):fz2(l)) + &
				FieldSources(l) % magnitude * GAUSS( t + FieldSources(l) % shift )
			
			CASE (3) ! Ez
				IF ( FieldSources(l) % update == 1 ) Ez(fx1(l):fx2(l),fy1(l):fy2(l),fz1(l):fz2(l)) = 0
				Ez(fx1(l):fx2(l),fy1(l):fy2(l),fz1(l):fz2(l)) = &
				Ez(fx1(l):fx2(l),fy1(l):fy2(l),fz1(l):fz2(l)) + &
				FieldSources(l) % magnitude * GAUSS( t + FieldSources(l) % shift )
		
		END SELECT
	
	END DO

END IF



!!! Excitation and readout !!!

IF (lumpedsource_flag) THEN

	IF (hardsource_flag) THEN
	
	    ! Hard source (Rs == 0)
	    SELECT CASE (SourceDir)
	        CASE (1) ! x dir
	            out(t) = Ex(SourceSubs(1),SourceSubs(2),SourceSubs(3))
	            Ex(SourceSubs(1),SourceSubs(2),SourceSubs(3)) = in(t)
	        CASE (2) ! y dir
	            out(t) = Ey(SourceSubs(1),SourceSubs(2),SourceSubs(3))
	            Ey(SourceSubs(1),SourceSubs(2),SourceSubs(3)) = in(t)
	        CASE (3) ! z dir
	            out(t) = Ez(SourceSubs(1),SourceSubs(2),SourceSubs(3))
	            Ez(SourceSubs(1),SourceSubs(2),SourceSubs(3)) = in(t)
	    END SELECT
	
	ELSE
	
	    ! Resistive source (Rs /= 0)
	    SELECT CASE (SourceDir)
	        CASE (1) ! x dir
	            Ex(SourceSubs(1),SourceSubs(2),SourceSubs(3)) = &
	            Ex(SourceSubs(1),SourceSubs(2),SourceSubs(3)) - in(t)
	            out(t) = Ex(SourceSubs(1),SourceSubs(2),SourceSubs(3))
	        CASE (2) ! y dir
	            Ey(SourceSubs(1),SourceSubs(2),SourceSubs(3)) = &
	            Ey(SourceSubs(1),SourceSubs(2),SourceSubs(3)) - in(t)
	            out(t) = Ey(SourceSubs(1),SourceSubs(2),SourceSubs(3))
	        CASE (3) ! z dir
	            Ez(SourceSubs(1),SourceSubs(2),SourceSubs(3)) = &
	            Ez(SourceSubs(1),SourceSubs(2),SourceSubs(3)) - in(t)
	            out(t) = Ez(SourceSubs(1),SourceSubs(2),SourceSubs(3))
	    END SELECT
	
	END IF

END IF



! DFT of NFFF surfaces
IF (nfff_flag) THEN

DO p = 1,freq_ff_size
! x-y plane
Ex_xy_dft(:,:,:,p) = Ex_xy_dft(:,:,:,p) &
                     + Wff(t,p)*Ex( x_nfff_1(1):x_nfff(2), &
                                    y_nfff_1(1):y_nfff_1(2), &
                                    z_nfff_1 )
Ey_xy_dft(:,:,:,p) = Ey_xy_dft(:,:,:,p) &
                     + Wff(t,p)*Ey( x_nfff_1(1):x_nfff_1(2), &
                                    y_nfff_1(1):y_nfff(2), &
                                    z_nfff_1 )
! y-z plane
Ey_yz_dft(:,:,:,p) = Ey_yz_dft(:,:,:,p) &
                     + Wff(t,p)*Ey( x_nfff_1, &
                                    y_nfff_1(1):y_nfff(2), &
                                    z_nfff_1(1):z_nfff_1(2) )
Ez_yz_dft(:,:,:,p) = Ez_yz_dft(:,:,:,p) &
                     + Wff(t,p)*Ez( x_nfff_1, &
                                    y_nfff_1(1):y_nfff_1(2), &
                                    z_nfff_1(1):z_nfff(2) )
! z-x plane
Ez_zx_dft(:,:,:,p) = Ez_zx_dft(:,:,:,p) &
                     + Wff(t,p)*Ez( x_nfff_1(1):x_nfff_1(2), &
                                    y_nfff_1, &
                                    z_nfff_1(1):z_nfff(2) )
Ex_zx_dft(:,:,:,p) = Ex_zx_dft(:,:,:,p) &
                     + Wff(t,p)*Ex( x_nfff_1(1):x_nfff(2), &
                                    y_nfff_1, &
                                    z_nfff_1(1):z_nfff_1(2) )
END DO

END IF


! Near field DFT
IF (field_flag) THEN

DO p = 1,freq_nf_size
DO k=1,zmax
DO j=1,ymax
DO i=1,xmax
    Ex_field(i,j,k,p) = CMPLX(  REAL( Ex_field(i,j,k,p) ) +  REAL( Wnf(t,p) )*Ex(i,j,k), &
                               AIMAG( Ex_field(i,j,k,p) ) + AIMAG( Wnf(t,p) )*Ex(i,j,k), DOUBLE )
    Ey_field(i,j,k,p) = CMPLX(  REAL( Ey_field(i,j,k,p) ) +  REAL( Wnf(t,p) )*Ey(i,j,k), &
                               AIMAG( Ey_field(i,j,k,p) ) + AIMAG( Wnf(t,p) )*Ey(i,j,k), DOUBLE )
    Ez_field(i,j,k,p) = CMPLX(  REAL( Ez_field(i,j,k,p) ) +  REAL( Wnf(t,p) )*Ez(i,j,k), &
                               AIMAG( Ez_field(i,j,k,p) ) + AIMAG( Wnf(t,p) )*Ez(i,j,k), DOUBLE )
END DO
END DO
END DO
END DO

END IF


!!! H component update !!!

! regular region
DO k=1,zmax_pml
DO j=1,ymax_pml
DO i=1,xmax_pml
    Hx(i,j,k) = Hx(i,j,k) + &
                Dbx(media(i,j,k))*( Ey(i,j,indz2(k)) - Ey(i,j,k) &
                                  - Ez(i,indy2(j),k) + Ez(i,j,k) )
    Hy(i,j,k) = Hy(i,j,k) + &
                Dby(media(i,j,k))*( Ez(indx2(i),j,k) - Ez(i,j,k) &
                                  - Ex(i,j,indz2(k)) + Ex(i,j,k) )
    Hz(i,j,k) = Hz(i,j,k) + &
                Dbz(media(i,j,k))*( Ex(i,indy2(j),k) - Ex(i,j,k) &
                                  - Ey(indx2(i),j,k) + Ey(i,j,k) )
END DO
END DO
END DO


! PML region
IF (pml_flag) THEN

DO k=1,zmax_pml
DO j=1,ymax_pml
DO i=1,pmax2
    Psi_Hy_x(i,j,k) = b_h(i)*Psi_Hy_x(i,j,k) &
        + c_h(i)*( Ez(indx_pml2(i),j,k) - Ez(indx_pml(i),j,k) )
    Hy(indx_pml(i),j,k) = Hy(indx_pml(i),j,k) &
        + Dby(media(indx_pml(i),j,k))*Psi_Hy_x(i,j,k)
    Psi_Hz_x(i,j,k) = b_h(i)*Psi_Hz_x(i,j,k) &
        + c_h(i)*( Ey(indx_pml2(i),j,k) - Ey(indx_pml(i),j,k) )
    Hz(indx_pml(i),j,k) = Hz(indx_pml(i),j,k) &
        - Dbz(media(indx_pml(i),j,k))*Psi_Hz_x(i,j,k)
END DO
END DO
END DO

DO k=1,zmax_pml
DO j=1,pmax2
DO i=1,xmax_pml
    Psi_Hz_y(i,j,k) = b_h(j)*Psi_Hz_y(i,j,k) &
        + c_h(j)*( Ex(i,indy_pml2(j),k) - Ex(i,indy_pml(j),k) )
    Hz(i,indy_pml(j),k) = Hz(i,indy_pml(j),k) &
        + Dbz(media(i,indy_pml(j),k))*Psi_Hz_y(i,j,k)
    Psi_Hx_y(i,j,k) = b_h(j)*Psi_Hx_y(i,j,k) &
        + c_h(j)*( Ez(i,indy_pml2(j),k) - Ez(i,indy_pml(j),k) )
    Hx(i,indy_pml(j),k) = Hx(i,indy_pml(j),k) &
        - Dbx(media(i,indy_pml(j),k))*Psi_Hx_y(i,j,k)
END DO
END DO
END DO

DO k=1,pmax2
DO j=1,ymax_pml
DO i=1,xmax_pml
    Psi_Hx_z(i,j,k) = b_h(k)*Psi_Hx_z(i,j,k) &
        + c_h(k)*( Ey(i,j,indz_pml2(k)) - Ey(i,j,indz_pml(k)) )
    Hx(i,j,indz_pml(k)) = Hx(i,j,indz_pml(k)) &
        + Dbx(media(i,j,indz_pml(k)))*Psi_Hx_z(i,j,k)
    Psi_Hy_z(i,j,k) = b_h(k)*Psi_Hy_z(i,j,k) &
        + c_h(k)*( Ex(i,j,indz_pml2(k)) - Ex(i,j,indz_pml(k)) )
    Hy(i,j,indz_pml(k)) = Hy(i,j,indz_pml(k)) &
        - Dby(media(i,j,indz_pml(k)))*Psi_Hy_z(i,j,k)
END DO
END DO
END DO

END IF



! Field sources
IF (fieldsource_flag) THEN

	DO l = 1,fieldsource_count
	
		SELECT CASE ( FieldSources(l) % component )
		
			CASE (4) ! Hx
				IF ( FieldSources(l) % update == 1 ) Hx(fx1(l):fx2(l),fy1(l):fy2(l),fz1(l):fz2(l)) = 0
				Hx(fx1(l):fx2(l),fy1(l):fy2(l),fz1(l):fz2(l)) = &
				Hx(fx1(l):fx2(l),fy1(l):fy2(l),fz1(l):fz2(l)) + &
				FieldSources(l) % magnitude * GAUSS( t + FieldSources(l) % shift )
			
			CASE (5) ! Hy
				IF ( FieldSources(l) % update == 1 ) Hy(fx1(l):fx2(l),fy1(l):fy2(l),fz1(l):fz2(l)) = 0
				Hy(fx1(l):fx2(l),fy1(l):fy2(l),fz1(l):fz2(l)) = &
				Hy(fx1(l):fx2(l),fy1(l):fy2(l),fz1(l):fz2(l)) + &
				FieldSources(l) % magnitude * GAUSS( t + FieldSources(l) % shift )
			
			CASE (6) ! Hz
				IF ( FieldSources(l) % update == 1 ) Hz(fx1(l):fx2(l),fy1(l):fy2(l),fz1(l):fz2(l)) = 0
				Hz(fx1(l):fx2(l),fy1(l):fy2(l),fz1(l):fz2(l)) = &
				Hz(fx1(l):fx2(l),fy1(l):fy2(l),fz1(l):fz2(l)) + &
				FieldSources(l) % magnitude * GAUSS( t + FieldSources(l) % shift )
		
		END SELECT
	
	END DO

END IF



! DFT of NFFF surfaces
IF (nfff_flag) THEN

DO p = 1,freq_ff_size
! x-y plane
Hy_xy_dft(:,:,:,p) = Hy_xy_dft(:,:,:,p) &
                      + Wff(t,p)*Hy( x_nfff_1(1):x_nfff(2), &
                                     y_nfff_1(1):y_nfff_1(2), &
                                     z_nfff )
Hx_xy_dft(:,:,:,p) = Hx_xy_dft(:,:,:,p) &
                      + Wff(t,p)*Hx( x_nfff_1(1):x_nfff_1(2), &
                                     y_nfff_1(1):y_nfff(2), &
                                     z_nfff )
Hy_xy_dft2(:,:,:,p) = Hy_xy_dft2(:,:,:,p) &
                      + Wff(t,p)*Hy( x_nfff_1(1):x_nfff(2), &
                                     y_nfff_1(1):y_nfff_1(2), &
                                     z_nfff_1 )
Hx_xy_dft2(:,:,:,p) = Hx_xy_dft2(:,:,:,p) &
                      + Wff(t,p)*Hx( x_nfff_1(1):x_nfff_1(2), &
                                     y_nfff_1(1):y_nfff(2), &
                                     z_nfff_1 )
! y-z plane
Hz_yz_dft(:,:,:,p) = Hz_yz_dft(:,:,:,p) &
                      + Wff(t,p)*Hz( x_nfff, &
                                     y_nfff_1(1):y_nfff(2), &
                                     z_nfff_1(1):z_nfff_1(2) )
Hy_yz_dft(:,:,:,p) = Hy_yz_dft(:,:,:,p) &
                      + Wff(t,p)*Hy( x_nfff, &
                                     y_nfff_1(1):y_nfff_1(2), &
                                     z_nfff_1(1):z_nfff(2) )
Hz_yz_dft2(:,:,:,p) = Hz_yz_dft2(:,:,:,p) &
                      + Wff(t,p)*Hz( x_nfff_1, &
                                     y_nfff_1(1):y_nfff(2), &
                                     z_nfff_1(1):z_nfff_1(2) )
Hy_yz_dft2(:,:,:,p) = Hy_yz_dft2(:,:,:,p) &
                      + Wff(t,p)*Hy( x_nfff_1, &
                                     y_nfff_1(1):y_nfff_1(2), &
                                     z_nfff_1(1):z_nfff(2) )
! z-x plane
Hx_zx_dft(:,:,:,p) = Hx_zx_dft(:,:,:,p) &
                      + Wff(t,p)*Hx( x_nfff_1(1):x_nfff_1(2), &
                                     y_nfff, &
                                     z_nfff_1(1):z_nfff(2) )
Hz_zx_dft(:,:,:,p) = Hz_zx_dft(:,:,:,p) &
                      + Wff(t,p)*Hz( x_nfff_1(1):x_nfff(2), &
                                     y_nfff, &
                                     z_nfff_1(1):z_nfff_1(2) )
Hx_zx_dft2(:,:,:,p) = Hx_zx_dft2(:,:,:,p) &
                      + Wff(t,p)*Hx( x_nfff_1(1):x_nfff_1(2), &
                                     y_nfff_1, &
                                     z_nfff_1(1):z_nfff(2) )
Hz_zx_dft2(:,:,:,p) = Hz_zx_dft2(:,:,:,p) &
                      + Wff(t,p)*Hz( x_nfff_1(1):x_nfff(2), &
                                     y_nfff_1, &
                                     z_nfff_1(1):z_nfff_1(2) )
END DO

END IF


! Probes
IF (probe_flag) THEN

    DO i = 1,probe_count
    
        SELECT CASE (probe_type(i))
            CASE (1)
                probes(t,i) = Ex(x_probe(i),y_probe(i),z_probe(i))
            CASE (2)
                probes(t,i) = Ey(x_probe(i),y_probe(i),z_probe(i))
            CASE (3)
                probes(t,i) = Ez(x_probe(i),y_probe(i),z_probe(i))
            CASE (4)
                probes(t,i) = Hx(x_probe(i),y_probe(i),z_probe(i))
            CASE (5)
                probes(t,i) = Hy(x_probe(i),y_probe(i),z_probe(i))
            CASE (6)
                probes(t,i) = Hz(x_probe(i),y_probe(i),z_probe(i))
        END SELECT
    
    END DO

END IF


! Energy computation
IF (energy_flag) THEN
IF (t == e_comp) THEN

    energy_now = 0
    DO k=1,zmax
    DO j=1,ymax
    DO i=1,xmax

        energy_now = energy_now + &
            Ex(i,j,k)**2/Cbx0(media(i,j,k)) + &
            Ey(i,j,k)**2/Cby0(media(i,j,k)) + &
            Ez(i,j,k)**2/Cbz0(media(i,j,k)) + &
            Hx(i,j,k)**2/Dbx(media(i,j,k)) + &
            Hy(i,j,k)**2/Dby(media(i,j,k)) + &
            Hz(i,j,k)**2/Dbz(media(i,j,k))

    END DO
    END DO
    END DO
    
    ! Excluding the feeding point
    IF (lumpedsource_flag) THEN
	    SELECT CASE (SourceDir)
	        CASE (1) ! x dir
	            energy_now = energy_now - &
	                Ex(SourceSubs(1),SourceSubs(2),SourceSubs(3))**2 &
	                /Cbx0(media(SourceSubs(1),SourceSubs(2),SourceSubs(3)))
	        CASE (2) ! y dir
	            energy_now = energy_now - &
	                Ey(SourceSubs(1),SourceSubs(2),SourceSubs(3))**2 &
	                /Cby0(media(SourceSubs(1),SourceSubs(2),SourceSubs(3)))
	        CASE (3) ! z dir
	            energy_now = energy_now - &
	                Ez(SourceSubs(1),SourceSubs(2),SourceSubs(3))**2 &
	                /Cbz0(media(SourceSubs(1),SourceSubs(2),SourceSubs(3)))
	    END SELECT
    END IF

    energy_now = energy_now*energy_mult
    energy(t) = energy_now
    e_comp = e_comp + e_step

    IF (energy_now > energy_max) THEN
        energy_max = energy_now
    END IF

    energy_percent = energy_now/energy_max*100

END IF
END IF


! Total energy termination condition
IF (drop_flag) drop_cond = (t > delay) .AND. (energy_percent < Drop_pct)


! Source energy termination condition
IF (drop2_flag) THEN

    IF (hardsource_flag) THEN
        source_now = in(t)**2
    ELSE
        source_now = out(t)**2
    END IF

    IF (source_now > source_max) source_max = source_now

    source_percent = source_now/source_max*100

    IF (source_percent < Drop2_pct) THEN
        source_count = source_count + 1
    ELSE
        source_count = 0
    END IF
    
    drop2_cond = (t > delay) .AND. (source_count > source_delay)
    
END IF


CALL CPU_TIME(time2)
overall = NINT((time2-time1)/t*Nt)

! displaying time step
PRINT '(A12, A, I0, A, I0, A, I0, A, I0, A, ES9.3, A, F0.2, A)', & 
    str_file, &
    '; n=', t, '/', Nt, &
    '; time=', (overall-NINT(time2-time1)), &
    '/', overall, 'sec; energy=', energy_now, &
    '/', energy_percent, '%'

! break on energy (total and source) below specified percent
IF (drop_flag .OR. drop2_flag) THEN
IF (drop_cond .AND. drop2_cond) THEN
    EXIT
END IF
END IF



END DO



DEALLOCATE( Cax, Cay, Caz, Cbx, Cby, Cbz, Dbx, Dby, Dbz, &
            media, in, indx1, indy1, indz1, indx2, indy2, indz2, &
            Ex, Ey, Ez, Hx, Hy, Hz, freq_nf )

IF (pml_flag) THEN
    DEALLOCATE( indx_pml, indy_pml, indz_pml, &
                indx_pml1, indy_pml1, indz_pml1, &
                indx_pml2, indy_pml2, indz_pml2, &
                b_e, b_h, c_e, c_h, &
                Psi_Ey_x, Psi_Ez_x, Psi_Ez_y, &
                Psi_Ex_y, Psi_Ex_z, Psi_Ey_z, &
                Psi_Hy_x, Psi_Hz_x, Psi_Hz_y, &
                Psi_Hx_y, Psi_Hx_z, Psi_Hy_z  )
END IF

IF (energy_flag) THEN
    DEALLOCATE( Cbx0, Cby0, Cbz0 )
END IF

IF (field_flag) THEN
    DEALLOCATE( Wnf )
END IF



CALL ELAPSED_TIME( partial_times(2) )



! Radiation pattern computation
IF (nfff_flag) THEN


! spatial geometric averaging of H-fields
WHERE ( Hx_xy_dft /= 0 ) Hx_xy_dft = Hx_xy_dft*SQRT( Hx_xy_dft2/Hx_xy_dft )
WHERE ( Hy_xy_dft /= 0 ) Hy_xy_dft = Hy_xy_dft*SQRT( Hy_xy_dft2/Hy_xy_dft )
WHERE ( Hy_yz_dft /= 0 ) Hy_yz_dft = Hy_yz_dft*SQRT( Hy_yz_dft2/Hy_yz_dft )
WHERE ( Hz_yz_dft /= 0 ) Hz_yz_dft = Hz_yz_dft*SQRT( Hz_yz_dft2/Hz_yz_dft )
WHERE ( Hz_zx_dft /= 0 ) Hz_zx_dft = Hz_zx_dft*SQRT( Hz_zx_dft2/Hz_zx_dft )
WHERE ( Hx_zx_dft /= 0 ) Hx_zx_dft = Hx_zx_dft*SQRT( Hx_zx_dft2/Hx_zx_dft )

DEALLOCATE( Hx_xy_dft2, Hy_xy_dft2, & 
            Hy_yz_dft2, Hz_yz_dft2, &
            Hz_zx_dft2, Hx_zx_dft2, Wff )


! Zeroing switched off planes

! x-y plane
IF (nfff_faces(5) == 0) THEN
    Ex_xy_dft(:,:,1,:) = 0
    Ey_xy_dft(:,:,1,:) = 0
    Hx_xy_dft(:,:,1,:) = 0
    Hy_xy_dft(:,:,1,:) = 0
END IF
IF (nfff_faces(6) == 0) THEN
    Ex_xy_dft(:,:,2,:) = 0
    Ey_xy_dft(:,:,2,:) = 0
    Hx_xy_dft(:,:,2,:) = 0
    Hy_xy_dft(:,:,2,:) = 0
END IF

! y-z plane
IF (nfff_faces(1) == 0) THEN
    Ey_yz_dft(1,:,:,:) = 0
    Ez_yz_dft(1,:,:,:) = 0
    Hy_yz_dft(1,:,:,:) = 0
    Hz_yz_dft(1,:,:,:) = 0
END IF
IF (nfff_faces(2) == 0) THEN
    Ey_yz_dft(2,:,:,:) = 0
    Ez_yz_dft(2,:,:,:) = 0
    Hy_yz_dft(2,:,:,:) = 0
    Hz_yz_dft(2,:,:,:) = 0
END IF

! z-x plane
IF (nfff_faces(3) == 0) THEN
    Ez_zx_dft(:,1,:,:) = 0
    Ex_zx_dft(:,1,:,:) = 0
    Hz_zx_dft(:,1,:,:) = 0
    Hx_zx_dft(:,1,:,:) = 0
END IF
IF (nfff_faces(4) == 0) THEN
    Ez_zx_dft(:,2,:,:) = 0
    Ex_zx_dft(:,2,:,:) = 0
    Hz_zx_dft(:,2,:,:) = 0
    Hx_zx_dft(:,2,:,:) = 0
END IF

            
ALLOCATE( Eth(resolution(1)+1,resolution(2),freq_ff_size), &
          Eph(resolution(1)+1,resolution(2),freq_ff_size), &
          x1_orig(nfff_size_x1), x2_orig(nfff_size_x2), &
          y1_orig(nfff_size_y1), y2_orig(nfff_size_y2), &
          z1_orig(nfff_size_z1), z2_orig(nfff_size_z2), &
          cos_x1(nfff_size_x1), cos_x2(nfff_size_x2), &
          cos_y1(nfff_size_y1), cos_y2(nfff_size_y2), &
          cos_z1(nfff_size_z1), cos_z2(nfff_size_z2), &
          sin_x1(nfff_size_x1), sin_x2(nfff_size_x2), &
          sin_y1(nfff_size_y1), sin_y2(nfff_size_y2), &
          sin_z1(nfff_size_z1), sin_z2(nfff_size_z2)  )
Eth = 0; Eph = 0;

! theta and phi steps
Dth =   pi/resolution(1)
Dph = 2*pi/resolution(2)

! wave number
alpha = 2*pi*freq_ff/c;

! half-timestep shift and eta multiplication
DO p = 1,freq_ff_size

    phase_factor = eta*EXP( CMPLX(0.0, pi*freq_ff(p)*Dt, DOUBLE) )

    Hx_xy_dft(:,:,:,p) = Hx_xy_dft(:,:,:,p)*phase_factor
    Hy_xy_dft(:,:,:,p) = Hy_xy_dft(:,:,:,p)*phase_factor
    Hy_yz_dft(:,:,:,p) = Hy_yz_dft(:,:,:,p)*phase_factor
    Hz_yz_dft(:,:,:,p) = Hz_yz_dft(:,:,:,p)*phase_factor
    Hz_zx_dft(:,:,:,p) = Hz_zx_dft(:,:,:,p)*phase_factor
    Hx_zx_dft(:,:,:,p) = Hx_zx_dft(:,:,:,p)*phase_factor

END DO

! surfaces of edge boundary elements are halved
! x-y plane
Ex_xy_dft(:,[1, nfff_size_y1],:,:) = Ex_xy_dft(:,[1, nfff_size_y1],:,:)/2
Ey_xy_dft([1, nfff_size_x1],:,:,:) = Ey_xy_dft([1, nfff_size_x1],:,:,:)/2
Hx_xy_dft([1, nfff_size_x1],:,:,:) = Hx_xy_dft([1, nfff_size_x1],:,:,:)/2
Hy_xy_dft(:,[1, nfff_size_y1],:,:) = Hy_xy_dft(:,[1, nfff_size_y1],:,:)/2
! y-z plane
Ey_yz_dft(:,:,[1, nfff_size_z1],:) = Ey_yz_dft(:,:,[1, nfff_size_z1],:)/2
Ez_yz_dft(:,[1, nfff_size_y1],:,:) = Ez_yz_dft(:,[1, nfff_size_y1],:,:)/2
Hy_yz_dft(:,[1, nfff_size_y1],:,:) = Hy_yz_dft(:,[1, nfff_size_y1],:,:)/2
Hz_yz_dft(:,:,[1, nfff_size_z1],:) = Hz_yz_dft(:,:,[1, nfff_size_z1],:)/2
! z-x plane
Ez_zx_dft([1, nfff_size_x1],:,:,:) = Ez_zx_dft([1, nfff_size_x1],:,:,:)/2
Ex_zx_dft(:,:,[1, nfff_size_z1],:) = Ex_zx_dft(:,:,[1, nfff_size_z1],:)/2
Hz_zx_dft(:,:,[1, nfff_size_z1],:) = Hz_zx_dft(:,:,[1, nfff_size_z1],:)/2
Hx_zx_dft([1, nfff_size_x1],:,:,:) = Hx_zx_dft([1, nfff_size_x1],:,:,:)/2

! geometrical center * 2 !!!
geom_center_2 = [ x_nfff(1) + x_nfff(2), &
                  y_nfff(1) + y_nfff(2), &
                  z_nfff(1) + z_nfff(2) ]

! difference vector between geometrical and phase centers
r_GP = geom_center_2*(Dx/2) - phase_center

! position vectors on surfaces - originals
x1_orig = ( [ (i, i=2*x_nfff(1)  ,2*x_nfff(2)  ,2) ] - geom_center_2(1) )*(Dx/2)
x2_orig = ( [ (i, i=2*x_nfff(1)+1,2*x_nfff(2)-1,2) ] - geom_center_2(1) )*(Dx/2)
y1_orig = ( [ (i, i=2*y_nfff(1)  ,2*y_nfff(2)  ,2) ] - geom_center_2(2) )*(Dx/2)
y2_orig = ( [ (i, i=2*y_nfff(1)+1,2*y_nfff(2)-1,2) ] - geom_center_2(2) )*(Dx/2)
z1_orig = ( [ (i, i=2*z_nfff(1)  ,2*z_nfff(2)  ,2) ] - geom_center_2(3) )*(Dx/2)
z2_orig = ( [ (i, i=2*z_nfff(1)+1,2*z_nfff(2)-1,2) ] - geom_center_2(3) )*(Dx/2)

! overall coefficient
coef = CMPLX(0.0, -alpha*(Dx**2)/(4*pi), DOUBLE)



DO n = 1,resolution(2)

    PRINT *, n

    ! phi angle
    ph = (n-1)*Dph
    sin_ph = SIN(ph)
    cos_ph = COS(ph)

    ! phi-vector components
    ph_x = -sin_ph
    ph_y = cos_ph
    ph_z = 0

    DO m = 1,resolution(1)+1

        ! theta angle
        th = (m-1)*Dth
        sin_th = SIN(th)
        cos_th = COS(th)

        ! r-vector components
        r_x  = sin_th*cos_ph
        r_y  = sin_th*sin_ph
        r_z  = cos_th

        ! theta-vector components
        th_x = cos_th*cos_ph
        th_y = cos_th*sin_ph
        th_z = -sin_th

        DO p = 1,freq_ff_size

            ! COS and SIN of position vectors
            cos_x1 = COS( (alpha(p)*r_x)*x1_orig )
            cos_x2 = COS( (alpha(p)*r_x)*x2_orig )
            cos_y1 = COS( (alpha(p)*r_y)*y1_orig )
            cos_y2 = COS( (alpha(p)*r_y)*y2_orig )
            cos_z1 = COS( (alpha(p)*r_z)*z1_orig )
            cos_z2 = COS( (alpha(p)*r_z)*z2_orig )

            sin_x1 = SIN( (alpha(p)*r_x)*x1_orig )
            sin_x2 = SIN( (alpha(p)*r_x)*x2_orig )
            sin_y1 = SIN( (alpha(p)*r_y)*y1_orig )
            sin_y2 = SIN( (alpha(p)*r_y)*y2_orig )
            sin_z1 = SIN( (alpha(p)*r_z)*z1_orig )
            sin_z2 = SIN( (alpha(p)*r_z)*z2_orig )


            ! x-y plane

            ! Hx, Ey
            temp_H = 0;
            temp_E = 0;
               jd =  nfff_size_y2
            DO j = 1,nfff_size_y2
               id =  nfff_size_x1
            DO i = 1,nfff_size_x1
                cos_temp1 = cos_x1(i)*cos_y2(j) - sin_x1(i)*sin_y2(j)
                sin_temp1 = sin_x1(i)*cos_y2(j) + cos_x1(i)*sin_y2(j)
                cos_temp2 = cos_temp1*cos_z1(1) - sin_temp1*sin_z1(1)
                sin_temp2 = cos_temp1*sin_z1(1) + sin_temp1*cos_z1(1)
                temp_H = temp_H + CMPLX( (  REAL(Hx_xy_dft(id,jd,2,p)) -  REAL(Hx_xy_dft(i,j,1,p)) )*cos_temp2  &
                                       + ( AIMAG(Hx_xy_dft(id,jd,2,p)) + AIMAG(Hx_xy_dft(i,j,1,p)) )*sin_temp2, & 
                                         ( AIMAG(Hx_xy_dft(id,jd,2,p)) - AIMAG(Hx_xy_dft(i,j,1,p)) )*cos_temp2  &
                                       - (  REAL(Hx_xy_dft(id,jd,2,p)) +  REAL(Hx_xy_dft(i,j,1,p)) )*sin_temp2, &
                                         DOUBLE )
                temp_E = temp_E + CMPLX( (  REAL(Ey_xy_dft(id,jd,2,p))  -  REAL(Ey_xy_dft(i,j,1,p))  )*cos_temp2  &
                                       + ( AIMAG(Ey_xy_dft(id,jd,2,p))  + AIMAG(Ey_xy_dft(i,j,1,p))  )*sin_temp2, &
                                         ( AIMAG(Ey_xy_dft(id,jd,2,p))  - AIMAG(Ey_xy_dft(i,j,1,p))  )*cos_temp2  &
                                       - (  REAL(Ey_xy_dft(id,jd,2,p))  +  REAL(Ey_xy_dft(i,j,1,p))  )*sin_temp2, &
                                         DOUBLE )
                id = id - 1
            END DO
            jd = jd - 1
            END DO

            Eth(m,n,p) = Eth(m,n,p) + th_y*temp_H + ph_x*temp_E
            Eph(m,n,p) = Eph(m,n,p) + ph_y*temp_H - th_x*temp_E

            ! Hy, Ex
            temp_H = 0;
            temp_E = 0;
               jd =  nfff_size_y1
            DO j = 1,nfff_size_y1
               id =  nfff_size_x2
            DO i = 1,nfff_size_x2
                cos_temp1 = cos_x2(i)*cos_y1(j) - sin_x2(i)*sin_y1(j)
                sin_temp1 = sin_x2(i)*cos_y1(j) + cos_x2(i)*sin_y1(j)
                cos_temp2 = cos_temp1*cos_z1(1) - sin_temp1*sin_z1(1)
                sin_temp2 = cos_temp1*sin_z1(1) + sin_temp1*cos_z1(1)
                temp_H = temp_H + CMPLX( (  REAL(Hy_xy_dft(id,jd,2,p)) -  REAL(Hy_xy_dft(i,j,1,p)) )*cos_temp2  &
                                       + ( AIMAG(Hy_xy_dft(id,jd,2,p)) + AIMAG(Hy_xy_dft(i,j,1,p)) )*sin_temp2, &
                                         ( AIMAG(Hy_xy_dft(id,jd,2,p)) - AIMAG(Hy_xy_dft(i,j,1,p)) )*cos_temp2  &
                                       - (  REAL(Hy_xy_dft(id,jd,2,p)) +  REAL(Hy_xy_dft(i,j,1,p)) )*sin_temp2, &
                                         DOUBLE )
                temp_E = temp_E + CMPLX( (  REAL(Ex_xy_dft(id,jd,2,p))  -  REAL(Ex_xy_dft(i,j,1,p))  )*cos_temp2  &
                                       + ( AIMAG(Ex_xy_dft(id,jd,2,p))  + AIMAG(Ex_xy_dft(i,j,1,p))  )*sin_temp2, &
                                         ( AIMAG(Ex_xy_dft(id,jd,2,p))  - AIMAG(Ex_xy_dft(i,j,1,p))  )*cos_temp2  &
                                       - (  REAL(Ex_xy_dft(id,jd,2,p))  +  REAL(Ex_xy_dft(i,j,1,p))  )*sin_temp2, &
                                         DOUBLE )
                id = id - 1
            END DO
            jd = jd - 1
            END DO

            Eth(m,n,p) = Eth(m,n,p) - th_x*temp_H - ph_y*temp_E
            Eph(m,n,p) = Eph(m,n,p) - ph_x*temp_H + th_y*temp_E


            ! y-z plane

            ! Hy, Ez
            temp_H = 0;
            temp_E = 0;
               kd =  nfff_size_z2
            DO k = 1,nfff_size_z2
               jd =  nfff_size_y1
            DO j = 1,nfff_size_y1
                cos_temp1 = cos_x1(1)*cos_y1(j) - sin_x1(1)*sin_y1(j)
                sin_temp1 = sin_x1(1)*cos_y1(j) + cos_x1(1)*sin_y1(j)
                cos_temp2 = cos_temp1*cos_z2(k) - sin_temp1*sin_z2(k)
                sin_temp2 = cos_temp1*sin_z2(k) + sin_temp1*cos_z2(k)
                temp_H = temp_H + CMPLX( (  REAL(Hy_yz_dft(2,jd,kd,p)) -  REAL(Hy_yz_dft(1,j,k,p)) )*cos_temp2  &
                                       + ( AIMAG(Hy_yz_dft(2,jd,kd,p)) + AIMAG(Hy_yz_dft(1,j,k,p)) )*sin_temp2, &
                                         ( AIMAG(Hy_yz_dft(2,jd,kd,p)) - AIMAG(Hy_yz_dft(1,j,k,p)) )*cos_temp2  &
                                       - (  REAL(Hy_yz_dft(2,jd,kd,p)) +  REAL(Hy_yz_dft(1,j,k,p)) )*sin_temp2, &
                                         DOUBLE )
                temp_E = temp_E + CMPLX( (  REAL(Ez_yz_dft(2,jd,kd,p))  -  REAL(Ez_yz_dft(1,j,k,p))  )*cos_temp2  &
                                       + ( AIMAG(Ez_yz_dft(2,jd,kd,p))  + AIMAG(Ez_yz_dft(1,j,k,p))  )*sin_temp2, &
                                         ( AIMAG(Ez_yz_dft(2,jd,kd,p))  - AIMAG(Ez_yz_dft(1,j,k,p))  )*cos_temp2  &
                                       - (  REAL(Ez_yz_dft(2,jd,kd,p))  +  REAL(Ez_yz_dft(1,j,k,p))  )*sin_temp2, &
                                         DOUBLE )
                jd = jd - 1
            END DO
            kd = kd - 1
            END DO

            Eth(m,n,p) = Eth(m,n,p) + th_z*temp_H + ph_y*temp_E
            Eph(m,n,p) = Eph(m,n,p) + ph_z*temp_H - th_y*temp_E

            ! Hz, Ey
            temp_H = 0;
            temp_E = 0;
               kd =  nfff_size_z1
            DO k = 1,nfff_size_z1
               jd =  nfff_size_y2
            DO j = 1,nfff_size_y2
                cos_temp1 = cos_x1(1)*cos_y2(j) - sin_x1(1)*sin_y2(j)
                sin_temp1 = sin_x1(1)*cos_y2(j) + cos_x1(1)*sin_y2(j)
                cos_temp2 = cos_temp1*cos_z1(k) - sin_temp1*sin_z1(k)
                sin_temp2 = cos_temp1*sin_z1(k) + sin_temp1*cos_z1(k)
                temp_H = temp_H + CMPLX( (  REAL(Hz_yz_dft(2,jd,kd,p)) -  REAL(Hz_yz_dft(1,j,k,p)) )*cos_temp2  &
                                       + ( AIMAG(Hz_yz_dft(2,jd,kd,p)) + AIMAG(Hz_yz_dft(1,j,k,p)) )*sin_temp2, &
                                         ( AIMAG(Hz_yz_dft(2,jd,kd,p)) - AIMAG(Hz_yz_dft(1,j,k,p)) )*cos_temp2  &
                                       - (  REAL(Hz_yz_dft(2,jd,kd,p)) +  REAL(Hz_yz_dft(1,j,k,p)) )*sin_temp2, &
                                         DOUBLE )
                temp_E = temp_E + CMPLX( (  REAL(Ey_yz_dft(2,jd,kd,p))  -  REAL(Ey_yz_dft(1,j,k,p))  )*cos_temp2  &
                                       + ( AIMAG(Ey_yz_dft(2,jd,kd,p))  + AIMAG(Ey_yz_dft(1,j,k,p))  )*sin_temp2, &
                                         ( AIMAG(Ey_yz_dft(2,jd,kd,p))  - AIMAG(Ey_yz_dft(1,j,k,p))  )*cos_temp2  &
                                       - (  REAL(Ey_yz_dft(2,jd,kd,p))  +  REAL(Ey_yz_dft(1,j,k,p))  )*sin_temp2, &
                                         DOUBLE )
                jd = jd - 1
            END DO
            kd = kd - 1
            END DO

            Eth(m,n,p) = Eth(m,n,p) - th_y*temp_H - ph_z*temp_E
            Eph(m,n,p) = Eph(m,n,p) - ph_y*temp_H + th_z*temp_E


            ! z-x plane

            ! Hz, Ex
            temp_H = 0;
            temp_E = 0;
               kd =  nfff_size_z1
            DO k = 1,nfff_size_z1
               id =  nfff_size_x2
            DO i = 1,nfff_size_x2
                cos_temp1 = cos_x2(i)*cos_y1(1) - sin_x2(i)*sin_y1(1)
                sin_temp1 = sin_x2(i)*cos_y1(1) + cos_x2(i)*sin_y1(1)
                cos_temp2 = cos_temp1*cos_z1(k) - sin_temp1*sin_z1(k)
                sin_temp2 = cos_temp1*sin_z1(k) + sin_temp1*cos_z1(k)
                temp_H = temp_H + CMPLX( (  REAL(Hz_zx_dft(id,2,kd,p)) -  REAL(Hz_zx_dft(i,1,k,p)) )*cos_temp2  &
                                       + ( AIMAG(Hz_zx_dft(id,2,kd,p)) + AIMAG(Hz_zx_dft(i,1,k,p)) )*sin_temp2, &
                                         ( AIMAG(Hz_zx_dft(id,2,kd,p)) - AIMAG(Hz_zx_dft(i,1,k,p)) )*cos_temp2  &
                                       - (  REAL(Hz_zx_dft(id,2,kd,p)) +  REAL(Hz_zx_dft(i,1,k,p)) )*sin_temp2, &
                                         DOUBLE )
                temp_E = temp_E + CMPLX( (  REAL(Ex_zx_dft(id,2,kd,p))  -  REAL(Ex_zx_dft(i,1,k,p))  )*cos_temp2  &
                                       + ( AIMAG(Ex_zx_dft(id,2,kd,p))  + AIMAG(Ex_zx_dft(i,1,k,p))  )*sin_temp2, &
                                         ( AIMAG(Ex_zx_dft(id,2,kd,p))  - AIMAG(Ex_zx_dft(i,1,k,p))  )*cos_temp2  &
                                       - (  REAL(Ex_zx_dft(id,2,kd,p))  +  REAL(Ex_zx_dft(i,1,k,p))  )*sin_temp2, &
                                         DOUBLE )
                id = id - 1
            END DO
            kd = kd - 1
            END DO

            Eth(m,n,p) = Eth(m,n,p) + th_x*temp_H + ph_z*temp_E
            Eph(m,n,p) = Eph(m,n,p) + ph_x*temp_H - th_z*temp_E

            ! Hx, Ez
            temp_H = 0;
            temp_E = 0;
               kd =  nfff_size_z2
            DO k = 1,nfff_size_z2
               id =  nfff_size_x1
            DO i = 1,nfff_size_x1
                cos_temp1 = cos_x1(i)*cos_y1(1) - sin_x1(i)*sin_y1(1)
                sin_temp1 = sin_x1(i)*cos_y1(1) + cos_x1(i)*sin_y1(1)
                cos_temp2 = cos_temp1*cos_z2(k) - sin_temp1*sin_z2(k)
                sin_temp2 = cos_temp1*sin_z2(k) + sin_temp1*cos_z2(k)
                temp_H = temp_H + CMPLX( (  REAL(Hx_zx_dft(id,2,kd,p)) -  REAL(Hx_zx_dft(i,1,k,p)) )*cos_temp2  &
                                       + ( AIMAG(Hx_zx_dft(id,2,kd,p)) + AIMAG(Hx_zx_dft(i,1,k,p)) )*sin_temp2, &
                                         ( AIMAG(Hx_zx_dft(id,2,kd,p)) - AIMAG(Hx_zx_dft(i,1,k,p)) )*cos_temp2  &
                                       - (  REAL(Hx_zx_dft(id,2,kd,p)) +  REAL(Hx_zx_dft(i,1,k,p)) )*sin_temp2, &
                                         DOUBLE )
                temp_E = temp_E + CMPLX( (  REAL(Ez_zx_dft(id,2,kd,p))  -  REAL(Ez_zx_dft(i,1,k,p))  )*cos_temp2  &
                                       + ( AIMAG(Ez_zx_dft(id,2,kd,p))  + AIMAG(Ez_zx_dft(i,1,k,p))  )*sin_temp2, &
                                         ( AIMAG(Ez_zx_dft(id,2,kd,p))  - AIMAG(Ez_zx_dft(i,1,k,p))  )*cos_temp2  &
                                       - (  REAL(Ez_zx_dft(id,2,kd,p))  +  REAL(Ez_zx_dft(i,1,k,p))  )*sin_temp2, &
                                         DOUBLE )
                id = id - 1
            END DO
            kd = kd - 1
            END DO

            Eth(m,n,p) = Eth(m,n,p) - th_z*temp_H - ph_x*temp_E
            Eph(m,n,p) = Eph(m,n,p) - ph_z*temp_H + th_x*temp_E


            ! phase of the G-P difference
            coef_GP = EXP( CMPLX(0.0, alpha(p)*( r_GP(1)*r_x + r_GP(2)*r_y + r_GP(3)*r_z ), DOUBLE) )

            ! final multiplication
            Eth(m,n,p) = Eth(m,n,p)*coef(p)*coef_GP
            Eph(m,n,p) = Eph(m,n,p)*coef(p)*coef_GP

        END DO ! end p

    END DO

END DO


END IF



CALL ELAPSED_TIME( partial_times(3) )



! Writing file
OPEN (22, FILE=(TRIM(str_file)//'.out'), STATUS='REPLACE', ACCESS='STREAM', RECL=HUGE(0_8))
WRITE (22) out, overall
IF (nfff_flag) WRITE (22) REAL(Eth), AIMAG(Eth), REAL(Eph), AIMAG(Eph)

! This is necessary so that the buffer limit (?) is not exceeded at write
IF (field_flag) THEN

    DO p = 1,freq_nf_size
    DO k = 1,zmax
        WRITE (22) REAL(Ex_field(1:xmax,1:ymax,k,p))
    END DO
    END DO

    DO p = 1,freq_nf_size
    DO k = 1,zmax
        WRITE (22) AIMAG(Ex_field(1:xmax,1:ymax,k,p))
    END DO
    END DO

    DO p = 1,freq_nf_size
    DO k = 1,zmax
        WRITE (22) REAL(Ey_field(1:xmax,1:ymax,k,p))
    END DO
    END DO

    DO p = 1,freq_nf_size
    DO k = 1,zmax
        WRITE (22) AIMAG(Ey_field(1:xmax,1:ymax,k,p))
    END DO
    END DO

    DO p = 1,freq_nf_size
    DO k = 1,zmax
        WRITE (22) REAL(Ez_field(1:xmax,1:ymax,k,p))
    END DO
    END DO

    DO p = 1,freq_nf_size
    DO k = 1,zmax
        WRITE (22) AIMAG(Ez_field(1:xmax,1:ymax,k,p))
    END DO
    END DO

END IF
        
IF (energy_flag) WRITE (22) energy

IF (probe_flag) WRITE (22) probes

WRITE (22) Rs

CLOSE (22)



! Deallocation of variables (necessary for cycling)
DEALLOCATE( out, freq_ff )
IF (nfff_flag) THEN
    DEALLOCATE( Ex_xy_dft, Ey_xy_dft, Hx_xy_dft, Hy_xy_dft, &
                Ey_yz_dft, Ez_yz_dft, Hy_yz_dft, Hz_yz_dft, &
                Ez_zx_dft, Ex_zx_dft, Hz_zx_dft, Hx_zx_dft, &
                alpha, coef, Eth, Eph, &
                x1_orig, x2_orig, y1_orig, y2_orig, z1_orig, z2_orig, &
                cos_x1, cos_x2, cos_y1, cos_y2, cos_z1, cos_z2, &
                sin_x1, sin_x2, sin_y1, sin_y2, sin_z1, sin_z2 )
END IF
IF (field_flag) THEN
    DEALLOCATE( Ex_field, Ey_field, Ez_field )
END IF
IF (energy_flag) THEN
    DEALLOCATE( energy )
END IF
IF (probe_flag) THEN
    DEALLOCATE( probes, x_probe, y_probe, z_probe, probe_type )
END IF
IF (fieldsource_flag) THEN
    DEALLOCATE( FieldSources, fx1, fx2, fy1, fy2, fz1, fz2 )
END IF



CALL DATE_AND_TIME( date, time, zone, time_end )
CALL ELAPSED_TIME( partial_times(4) )
    
! writing the times at the end of the file
OPEN (22, FILE=(TRIM(str_file)//'.out'), STATUS='OLD', ACCESS='STREAM', FORM='UNFORMATTED', RECL=HUGE(0_8))
INQUIRE (22, SIZE=filesize)
WRITE (22, POS=filesize+1) time_start, time_end, partial_times
CLOSE (22)



END DO





CONTAINS

    SUBROUTINE ELAPSED_TIME( time_arg )
    
        IMPLICIT NONE
    
        INTEGER, INTENT(OUT) :: time_arg
        INTEGER :: count_clock, count_rate, new_time, last_time = 0
        
        CALL SYSTEM_CLOCK( count_clock, count_rate )
        
        IF (count_rate == 0) THEN
            time_arg = 0
        ELSE
            new_time = count_clock / count_rate
            time_arg = new_time - last_time
            last_time = new_time
        END IF
    
    END SUBROUTINE ELAPSED_TIME
    
    
    
    SUBROUTINE STAIRCASEWIRE( point1, point2 )
    
        IMPLICIT NONE

        REAL(DOUBLE), INTENT(IN) :: point1(3), point2(3)
    
        INTEGER(4) :: steps(3), dir(3), cr(3), act_step(3), ind, ind_temp(1), &
                      round_point1(3), round_point2(3), act_point(3)
        REAL(DOUBLE) :: vector(3), offset(3), span(3), parameters(3)
        
        vector = point2 - point1
        
        round_point1 = NINT( point1, 4 )
        round_point2 = NINT( point2, 4 )
        
        steps = ABS( round_point2 - round_point1 )
        
        cr = SIGN( [1,2,3], NINT( vector ) ) ! cr is negative where opposite direction
        WHERE ( steps == 0 ) cr = ABS(cr) ! put it positive in direction where there are no steps
        dir = SIGN( 1, cr ) ! has -1 where opposite direction, 1 elsewhere

        ! Simple bubble sort
        IF ( cr(1) > cr(2) ) cr = cr([2,1,3])
        IF ( cr(2) > cr(3) ) cr = cr([1,3,2])
        IF ( cr(1) > cr(2) ) cr = cr([2,1,3])
        
        cr = ABS(cr) ! now the directions are sorted with the opposite first [-3,-2,-1,1,2,3]
        
        offset = dir*( round_point1 - point1 ) ! offset of the first point from the grid
        span = ABS( vector ) 
        
        subs_length = SUM( steps )
        ALLOCATE( subs( subs_length, 3 ), comp( subs_length ) )
        
        act_step = 1
        
        WHERE ( span == 0 )
            parameters = HUGE( 0_DOUBLE )
        ELSEWHERE
            parameters = ( act_step - 0.5_DOUBLE + offset )/span
        END WHERE
        
        act_point = round_point1 - limits_1 + 1 ! normalized to media domain
        
        DO i = 1,subs_length
        
            ind_temp = cr( MINLOC( parameters(cr) ) )
            ind = ind_temp(1)
            subs(i,:) = act_point
            IF ( dir(ind) < 0 ) subs(i,ind) = subs(i,ind) - 1
            comp(i) = ind * dir(ind)
            act_point(ind) = act_point(ind) + dir(ind)
            act_step(ind)  = act_step(ind)  + 1
            parameters(ind) = ( act_step(ind) - 0.5_DOUBLE + offset(ind) )/span(ind)
        
        END DO
                    
    END SUBROUTINE STAIRCASEWIRE
    
    
    
    FUNCTION CROSS( vector1, vector2 )
    
        IMPLICIT NONE

        REAL(DOUBLE) :: cross(3)
        REAL(DOUBLE), INTENT(IN)  :: vector1(3), vector2(3)
        
        cross(1) = vector1(2) * vector2(3) - vector1(3) * vector2(2)
        cross(2) = vector1(3) * vector2(1) - vector1(1) * vector2(3)
        cross(3) = vector1(1) * vector2(2) - vector1(2) * vector2(1)
    
    END FUNCTION CROSS
    
    
    
    SUBROUTINE SEARCH_LUT
    
        IMPLICIT NONE

        mat_number = materials + 1 ! preliminary new material
        
        ! check if among the old ones
        ! check only this wire and in reverse order - higher likelihood of a hit
        DO n = materials,materials_prev+1,-1
            IF ( lut_actual(1) == lut(n,1) .AND. &
                 lut_actual(2) == lut(n,2) .AND. &
                 lut_actual(3) == lut(n,3) .AND. &
                 lut_actual(4) == lut(n,4) .AND. &
                 lut_actual(5) == lut(n,5) .AND. &
                 lut_actual(6) == lut(n,6) .AND. &
                 lut_actual(7) == lut(n,7) .AND. &
                 lut_actual(8) == lut(n,8) .AND. &
                 lut_actual(9) == lut(n,9) ) THEN
                mat_number = n
                EXIT
            END IF
        END DO
        
        IF ( mat_number > materials ) THEN
        
            materials = materials + 1
            IF ( materials > max_length ) STOP 'Error -- lookup table length exceeded'
            lut(mat_number,:) = lut_actual
            
        END IF

    END SUBROUTINE SEARCH_LUT

    
    
    SUBROUTINE ADD_TWEAKS( i )
    
        IMPLICIT NONE

        INTEGER(4), INTENT(IN) :: i
        
        DO j = 1,tweaks_sizes(i)
        
            READ (11) point2, tweak1, tweak2 
            round_point2 = NINT( point2, 4 ) + round_point1
            
            IF ( ANY( round_point2 < 1 ) .OR. ANY( round_point2 > [xmax,ymax,zmax] ) ) CYCLE
            
            mat_number = media(round_point2(1),round_point2(2),round_point2(3))
            lut_actual = lut(mat_number,:)
            IF ( i < 4 .AND. lut_actual(i) == 0 ) CYCLE ! the eps zero values are of PEC - fixed
            lut_actual(i) = lut_actual(i) * tweak2 + tweak1
            
            IF ( mat_number > materials_prev ) THEN
                lut(mat_number,:) = lut_actual
            ELSE
                materials = materials + 1
                IF ( materials > max_length ) STOP 'Error -- lookup table length exceeded'
                lut(materials,:) = lut_actual
                media(round_point2(1),round_point2(2),round_point2(3)) = materials
            END IF
        
        END DO
    
    END SUBROUTINE ADD_TWEAKS
    
    
    ELEMENTAL FUNCTION GAUSS( timestep )
    
        IMPLICIT NONE

        REAL(DOUBLE) :: gauss, timestep_shifted
        REAL(DOUBLE), INTENT(IN) :: timestep

        timestep_shifted = timestep - 1 - offset_dt
        
        IF ( ( timestep_shifted < -offset_dt ) .OR. ( timestep_shifted > offset_dt ) ) THEN
        	gauss = 0
        ELSE
        	gauss = mag * SIN( (2*pi*fc*Dt) * timestep_shifted ) &
        				* EXP( -(( (fs*Dt) * timestep_shifted )**2) )
        END IF
    
    END FUNCTION GAUSS


END PROGRAM FDTD



SUBROUTINE INTERRUPT

USE SHARED

IMPLICIT NONE

CHARACTER :: answer

DO
    WRITE (*, '(A)', ADVANCE='NO') 'Program interrupted. Abort/Continue? >'
    READ (*, '(A)') answer
    IF (answer=='c') RETURN
    IF (answer=='a') EXIT
END DO

DO
    WRITE (*, '(A)', ADVANCE='NO') 'Save unfinished results? Yes/No? >'
    READ (*, '(A)') answer
    IF (answer=='n') STOP
    IF (answer=='y') EXIT
END DO

WRITE (*, '(A)', ADVANCE='NO') 'Saving... '

! Writing file
OPEN (22, FILE='abort.out', STATUS='REPLACE', ACCESS='STREAM', RECL=HUGE(0_8))
WRITE (22) out, overall
IF (nfff_flag) THEN
    IF (.NOT.ALLOCATED(Eth)) THEN
        ALLOCATE(Eth(resolution(1)+1,resolution(2),freq_ff_size))
        Eth = 0;
    END IF
    IF (.NOT.ALLOCATED(Eph)) THEN
        ALLOCATE(Eph(resolution(1)+1,resolution(2),freq_ff_size))
        Eph = 0;
    END IF
    WRITE (22) REAL(Eth), AIMAG(Eth), REAL(Eph), AIMAG(Eph)
END IF
IF (field_flag) WRITE (22) REAL(Ex_field(1:xmax,1:ymax,1:zmax,:)), AIMAG(Ex_field(1:xmax,1:ymax,1:zmax,:)), &
                      REAL(Ey_field(1:xmax,1:ymax,1:zmax,:)), AIMAG(Ey_field(1:xmax,1:ymax,1:zmax,:)), &
                      REAL(Ez_field(1:xmax,1:ymax,1:zmax,:)), AIMAG(Ez_field(1:xmax,1:ymax,1:zmax,:))
IF (energy_flag) WRITE (22) energy
CLOSE (22)

WRITE (*, '(A)') 'done.'

STOP

END SUBROUTINE INTERRUPT
      
