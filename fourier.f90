!234567
       PROGRAM fourier

       implicit none
       integer, Parameter :: rangt = 10000, rangw = 10000
       Real*8 f(rangt), t(rangt), reftransform, imftransform
       Real*8 real_fourier(rangw), im_fourier(rangw), w(rangw)
       Real*8 real_function(rangt), im_function(rangt)
       Real*8 anti_reftransform, anti_imftransform
       Real*8, Parameter :: ti = -2., tf = 2., pi = 3.14
       Real*8, Parameter :: wi = -6., wf = 6.
       integer i, j, k

       open(unit = 2, file = 'fourier.dat')
       open(unit = 3, file = 'fouriert.dat')
       open(unit = 4, file = 'anti_fourier.dat')

       do i = 1, rangt
         f(i) = 0.0
         t(i) = 0.0
       enddo
       
       do i = 1, rangw
         real_fourier(i) = 0.0
         im_fourier  (i) = 0.0
         w(i) = wi + (wf-wi)/(rangw-1)*real(i-1)     
       enddo

       do i = 1, rangt
         t(i) = ti + real(i-1)*(tf-ti)/(rangt-1) 
       enddo

       do i = 1, rangt
         f(i) = cos(2.0*pi*(3.0*t(i)))*exp(-pi*t(i)**2)
       enddo

       do i = 1, rangt
         write(2,*) t(i), f(i)
       enddo
       
       do j = 1, rangw
         call ftransform(f,t,rangt,rangw,w,j,reftransform,imftransform)
         real_fourier(j) = reftransform
         im_fourier  (j) = imftransform
       enddo

!       do j = 1, rangt
!         call anti_ftransform()

       do j = 1, rangw
         write(3,*) w(j), real_fourier(j), im_fourier(j) 
       enddo

       do j = 1, rangt
         call  anti_ftransform(real_fourier,im_fourier,rangw,rangt&
                                   ,w,t,j,                            &
                                   anti_reftransform,anti_imftransform)
         real_function(j) = anti_reftransform
         im_function  (j) = anti_imftransform
       enddo

       do j = 1, rangt
         write(4,*) t(j), real_function(j), im_function(j)
       enddo

       end program



!======================================================================
!=                         FOURIER TRANSFORM                          =
!======================================================================

       subroutine  ftransform(f,t,rangt,rangw,w,j,&
                              reftransform,imftransform)
       implicit none
       
       Real*8, Parameter :: pi = 3.14
       Real*8 f(rangt), t(rangt), dt, reftransform, imftransform
       Real*8 w(rangw)
       integer rangt, rangw, i, j

       reftransform = 0.0
       dt         = (t(rangt)-t(1))/(rangt-1)
       
       do i = 2, rangt-1
         reftransform = reftransform + (f(i)*cos(2.*pi*w(j)*t(i)))*dt
       enddo
        
       reftransform = reftransform + (f(1)*cos(2.*pi*w(j)*t(1))/2.0 + &
       f(rangt)*cos(2.*pi*w(j)*t(rangt))/2.0)*dt



       imftransform = 0.0

       do i = 2, rangt-1
         imftransform = imftransform + (f(i)*sin(-2.*pi*w(j)*t(i)))*dt
       enddo

       imftransform = imftransform + (f(1)*sin(-2.*pi*w(j)*t(1))/2.0+ &
       f(rangt)*sin(-2.*pi*w(j)*t(rangt))/2.0)*dt
 
       return
       end subroutine ftransform


!======================================================================
!=                         FOURIER ANTITRANSFORM                      =
!======================================================================

       subroutine  anti_ftransform(real_fourier,im_fourier,rangw,rangt&
                                   ,w,t,j,                            &
                                   anti_reftransform,anti_imftransform)
       implicit none
       
       Real*8, Parameter :: pi = 3.14
       Real*8 w(rangw), dw, anti_reftransform, anti_imftransform
       Real*8 real_fourier(rangw), im_fourier(rangw), t(rangt)
       integer rangw, rangt, i, j

       anti_reftransform = 0.0
       dw         = (w(rangw)-w(1))/(rangw-1)
       
       do i = 2, rangw-1
         anti_reftransform = anti_reftransform + &
                           (real_fourier(i)*cos(2.*pi*w(i)*t(j)))*dw +&
                           (-im_fourier(i)*sin(2.*pi*w(i)*t(j)))*dw
       enddo
        
       anti_reftransform = anti_reftransform + &
                 (real_fourier(1)*cos(2.*pi*w(1)*t(j))/2.0 + &
                 real_fourier(rangw)*cos(2.*pi*w(rangw)*t(j))/2.0)*dw+&
                 (-im_fourier(1)*sin(2.*pi*w(1)*t(j))/2.0 - &
                 im_fourier(rangw)*sin(2.*pi*w(rangw)*t(j))/2.0)*dw




       anti_imftransform = 0.0

       do i = 2, rangw-1
         anti_imftransform = anti_imftransform + &
                           (real_fourier(i)*sin(2.*pi*w(i)*t(j)))*dw +&
                           (im_fourier(i)*cos(2.*pi*w(i)*t(j)))*dw
       enddo

       anti_imftransform = anti_imftransform + &
                 (real_fourier(1)*sin(2.*pi*w(1)*t(j))/2.0 + &
                 real_fourier(rangw)*sin(2.*pi*w(rangw)*t(j))/2.0)*dw+&
                 (im_fourier(1)*cos(2.*pi*w(1)*t(j))/2.0 + &
                 im_fourier(rangw)*cos(2.*pi*w(rangw)*t(j))/2.0)*dw


 
       return
       end subroutine anti_ftransform







