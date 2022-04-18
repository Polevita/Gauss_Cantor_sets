module dimsetup_subs
    integer, parameter :: short = SELECTED_INT_KIND(2)
    integer, parameter :: long = SELECTED_INT_KIND(12)
    integer, parameter :: qp = selected_real_kind(33, 4931)


contains
    
    function tempmatrix(temp, alphbet, n) result(mm) 
        implicit none 
        integer, intent(in) :: n 
        integer, intent(in) :: temp(:), alphbet(:) 

        integer(long) :: mm(4), xx(4)

        integer k0, k1, k2 

        mm(1) = 1
        mm(2) = 0 
        mm(3) = 0 
        mm(4) = 1 

        do k0 = 1, n
            xx = mm 
            mm(1) = xx(3)
            mm(2) = xx(4) 
            mm(3) = alphbet(temp(n-k0+1))*xx(3)+xx(1)
            mm(4) = alphbet(temp(n-k0+1))*xx(4)+xx(2)
        enddo 

    endfunction tempmatrix



    function firstsearch(badwords, ex, N, alphbet, An, wordsfile, matfile) result(total)
        implicit none 
        integer, intent(in) :: badwords(:,:), alphbet(:), An, ex, N 
        character(100), intent(in) :: wordsfile, matfile
        
        integer :: total 
        integer :: j, k, k0, k1, k2, k3
        integer, allocatable :: temp(:) 
        integer(long) :: mm(4)
        character(25) :: myprint, mymatrix
        logical :: sw0, sw1, sw2
        real(qp) :: x0, x1, x2, x3

        real*4 sec, sec1

        x0 = real(maxval(alphbet),qp)
!        print *, x0, 'x0'
        x1 = (x0 + sqrt(x0*x0 + 4.0_qp))*0.50_qp
!        print *, x1, 'x1'
        x2 = 0.50_qp*(x0 + (x0*x0 + 2.0_qp)/(sqrt(x0*x0+4.0_qp))) 
!        print *, x2, 'x2'
       
        k0 = ceiling(N*log10(x1)+log(x2))+1
!        print *, k0, 'maximal matrix entry'

        if (RANGE(mm(2))<k0) then
            print *, k0, 'k0', RANGE(mm(2)), ' -- Matrix. Error. Not enough space for the matrix. Stopping'
            stop
        endif 

        allocate(temp(N))
        sec = secnds(0.0)
        open(7, file = trim(wordsfile), status = 'replace')
        write (myprint, '(A1I4A7)') trim("("), N, trim("(1xI4))")

        open(9, file = trim(matfile), status = 'replace')
        write (mymatrix, '(A6I2A2)') trim("(4(4xI"), k0, trim("))")
        !print *, mymatrix

        total = 0 
        mainloop: do j = 0, An**N-1
            temp = 0
            k = j
            modloop: do k1 = 1,N 
                if (k>=1) then 
                    temp(N+1-k1) = mod(k,An) 
                    k = floor(real(k)/real(An))
                else 
                    exit modloop
                endif
            enddo modloop
            temp = temp+1 
            sw0 = .false.
            searchloop: do k = 1,ex 
                sw1 = .false.
                if (badwords(k,1) < N+1)  then 
                    wordloop: do k2 = 1, N + 1 - badwords(k,1)
                        sw2 = .false. 
                        checkloop: do k3 = 0, badwords(k,1) -1
                            if (temp(k2+k3) .ne. badwords(k,k3+2) ) then 
                                sw2 = .true.  ! good words
                                exit checkloop
                            endif
                        enddo checkloop
                        if (.not.sw2) then
                            sw1 = .true.  !bad word 
                            exit wordloop
                        endif
                    enddo wordloop
                endif
                if (sw1) then
                    sw0 = .true. !bad word
                    exit searchloop
                endif  
            enddo searchloop
            if (.not.sw0) then 
                write (7, FMT = myprint)  (temp(k), k=1, N )    
                mm = tempmatrix(temp,alphbet,N)
                write (9, FMT = mymatrix) (mm(k), k = 1, 4)    
                total = total + 1
            endif

        enddo mainloop

        sec1 = secnds(sec)
        print *, 'Words in total: ', total, '; : ', sec1, 'seconds to run'

        close (7)
        close (9)

        deallocate(temp)

    end function firstsearch

    subroutine wordsearch(total,N,ex,badwords,wordsfile,mulfile,comulfile,allwords,maxmul,uni,couni)
        implicit none 

        integer, intent(in) :: badwords(:,:) 
        character(100), intent(in) :: wordsfile, mulfile, comulfile
        integer, intent(in) :: total, N, ex
        integer, intent(out) :: maxmul, uni, couni, allwords(:,:) 

        integer, allocatable :: swn(:,:),  temp(:)
        logical, allocatable :: swords(:,:), indx(:) 
        logical :: sw0, sw1, sw2
        integer :: j, k, k0, k1, k2, k3, num, wmax, cmax

        character(20) :: myprint

        real*4 sec, sec1

        wmax = total
        cmax = total

        open(12,file = trim(wordsfile))

        do j = 1, total 
            read(12, *) (allwords(j,k), k = 1,N)
        enddo

        allocate(indx(total))
        allocate(temp(2*N))
        allocate(swords(wmax,total))
        allocate(swn(wmax,cmax))  
        sec = secnds(0.0)

        num = 0
        swn = 0 
        temp = 0
        mainloop: do j = 1, total
            indx = .false. 
            temp(N+1:2*N) = allwords(j,:) 
            do k = 1,  total 
                temp(1:N) = allwords(k,:)
                sw0 = .false.
                seqloop: do k0 = 1,ex
                    sw1 = .false.
                    mwloop: do k1 = N-badwords(k0,1)+2,N
                        sw2 = .false.
                        wordloop: do k2 = 0,badwords(k0,1)-1
                            if (temp(k1+k2).ne.badwords(k0,k2+2)) then
                                sw2 = .true. ! allowed
                                !print *, 'k0 =', k0, 'YES', (temp(k3), k3 = 1,N), 'AND', (temp(k3), k3 = N+1, 2*N)
                                exit wordloop
                            endif
                        enddo wordloop
                        if (sw2.eqv..false.) then
                            sw1 = .true. ! forbidden
                            !print *, 'k0 =', k0, 'NAY', (temp(k3), k3 = 1,N), 'AND', (temp(k3), k3 = N+1, 2*N)
                            exit mwloop
                        endif
                    enddo mwloop
                    if (sw1.eqv..true.) then
                        sw0 = .true. ! forbidden 
                        exit seqloop
                    endif
                enddo seqloop
                indx(k) = .not.sw0
            enddo
            if (j == 1) then 
                num = num+1
                swords(num,1:total) = indx
                swn(num,2) = 1;
                swn(num,1) = 3;
            else
                sw0 = .false.
                checkloop: do k = 1,num
                    if ((ALL(indx.eqv.swords(k,:))).eqv..true.) then 
                        sw0 = .true. 
                        swn(k,swn(k,1)) = j
                        swn(k,1) = swn (k,1) + 1    
                        exit checkloop
                    endif
                enddo checkloop
                if (sw0.eqv..false.) then 
                    num = num+1
                    swords(num,:) = indx
                    swn(num,2) = j
                    swn(num,1) = 3
                endif
            endif

            if (num > wmax) then 
                print *, 'too many words - something wrong'
                exit mainloop
            endif

            if (maxval(swn(:,1)) > cmax ) then 
                print *, 'too many coincidences'
                exit mainloop
            endif  

            if (mod(j,1000) == 0) then 
                print *, 'parsed', j, 'words'
            endif
        enddo mainloop

        open(7, file = trim(mulfile), status = 'replace')

        do j = 1, num
            write (myprint, '(A1I6A7)') trim("("), swn(j,1)-1, trim("(1xI6))")
            write (7, FMT = myprint)  (swn(j,k), k=1, (swn(j,1)-1) )    
        enddo

        close (7)

        uni = num;
        maxmul = maxval(swn(:,1))

        num = 0
        swn = 0 
        temp = 0
        swords = .false.

        mainloop2: do j = 1, total
            indx = .false. 
            temp(1:N) = allwords(j,:) 
            do k = 1,  total 
                temp(N+1:2*N) = allwords(k,:)
                sw0 = .false.
                seqloop2: do k0 = 1,ex
                    sw1 = .false.
                    mwloop2: do k1 = N-badwords(k0,1)+2,N
                        sw2 = .false.
                        wordloop2: do k2 = 0,badwords(k0,1)-1
                            if (temp(k1+k2).ne.badwords(k0,k2+2)) then
                                sw2 = .true. ! allowed
                                ! print *, 'k0 =', k0, 'YES', (temp(k3), k3 = 1,N), 'AND', (temp(k3), k3 = N+1, 2*N)
                                exit wordloop2
                            endif
                        enddo wordloop2
                        if (sw2.eqv..false.) then
                            sw1 = .true. ! forbidden
                            ! print *, 'k0 =', k0, 'NAY', (temp(k3), k3 = 1,N), 'AND', (temp(k3), k3 = N+1, 2*N)
                            exit mwloop2
                        endif
                    enddo mwloop2
                    if (sw1.eqv..true.) then
                        sw0 = .true. ! forbidden 
                        exit seqloop2
                    endif
                enddo seqloop2
                indx(k) = .not.sw0
            enddo
            if (j == 1) then 
                num = num+1
                swords(num,1:total) = indx
                swn(num,2) = 1;
                swn(num,1) = 3;
            else
                sw0 = .false.
                checkloop2: do k = 1,num
                    if ((ALL(indx.eqv.swords(k,:))).eqv..true.) then 
                        sw0 = .true. 
                        swn(k,swn(k,1)) = j
                        swn(k,1) = swn (k,1) + 1    
                        exit checkloop2
                    endif
                enddo checkloop2
                if (sw0.eqv..false.) then 
                    num = num+1
                    swords(num,:) = indx
                    swn(num,2) = j
                    swn(num,1) = 3
                endif
            endif

            if (num > wmax) then 
                print *, 'too many words - something wrong'
                exit mainloop2
            endif

            if (maxval(swn(:,1)) > cmax ) then 
                print *, 'too many coincidences'
                exit mainloop2
            endif  

            if (mod(j,1000) == 0) then 
                print *, 'parsed', j, 'words'
            endif
        enddo mainloop2

        open(7, file = trim(comulfile), status = 'replace')

        do j = 1, num
            write (myprint, '(A1I6A7)') trim("("), swn(j,1)-1, trim("(1xI6))")
            write (7, FMT = myprint)  (swn(j,k), k=1, (swn(j,1)-1) )    
        enddo

        close (7)
        
        couni = num 
        maxmul = max(maxmul, maxval(swn(:,1)))

        sec1 = secnds(sec)
        print *, 'Unique rows:', j-1, ' : ', sec1, 'seconds to run'

        deallocate(indx)
        deallocate(temp)
        deallocate(swords)
        deallocate(swn)  
    end subroutine wordsearch

    subroutine minimatrix(total,N,ex,maxmul,uni,couni,allwords,badwords,multi,comulti,outfile)
    implicit none 

    integer, intent(in) :: badwords(:,:), allwords(:,:), total, N, uni, couni, ex, maxmul
    character(100), intent(in) ::  multi, comulti, outfile !, muldata, comuldata
    
    integer(short), allocatable :: indx(:) 

    integer, allocatable :: uind(:), couind(:), temp(:), mul(:), comul(:)  
    logical :: sw0, sw1, sw2

    integer :: j, k, k0, k1, k2, k3
    character(25) :: myprint

    real*4 sec, sec1

    allocate(uind(uni))
    allocate(couind(couni)) 
    allocate(temp(maxmul))
    allocate(mul(total))
    allocate(comul(total))

    open(14, file = trim(comulti))
    open(12, file = trim(multi))

    mul = 0
    comul = 0
    do j = 1, uni 
        temp = 0
        read(12, *) temp(1), (temp(k), k = 2, temp(1)-1)
        uind(j) = temp(2)
        do k = 2, temp(1)-1
            mul(temp(k)) = j
        enddo
    enddo
    do j = 1, couni 
        temp = 0
        read(14, *) temp(1), (temp(k), k = 2, temp(1)-1)
        couind(j) = temp(2)
        do k = 2, temp(1)-1
            comul(temp(k)) = j
        enddo
    enddo
    deallocate(temp)
    close(12)
    close(14)

    if (minval(mul) < 1) then 
        print *, 'error: zero element in multiplicativity list. stopping.'
        stop
    endif   

!    write (myprint, '(A1I2A8)') trim("("), total, trim("(I4 /))")
!    print *, myprint
!    open(7, file = trim(muldata), status = 'replace')
!        write (7, FMT = myprint)  (mul(k), k=1,total)    
!    close(7)
    
!    open(9, file = trim(comuldata), status = 'replace')
!        write (9, FMT = myprint)  (comul(k), k=1,total)    
!    close(9)

    deallocate(mul)
    deallocate(comul)

    allocate(indx(couni))
    allocate(temp(2*N))

    sec = secnds(0.0)
    open(7, file = trim(outfile), status = 'replace')
    write (myprint, '(A1I2A7)') trim("("), uni, trim("(1xI3))")
    !print *, myprint

    mainloop: do j = 1, uni
        indx = 0 !.false. 
        temp(N+1:2*N) = allwords(uind(j),:) 
        do k = 1, couni 
            temp(1:N) = allwords(couind(k),:)
            sw0 = .false.
            seqloop: do k0 = 1,ex
                sw1 = .false.
                mwloop: do k1 = N-badwords(k0,1)+2,N
                    sw2 = .false.
                    wordloop: do k2 = 0,badwords(k0,1)-1
                        if (temp(k1+k2).ne.badwords(k0,k2+2)) then
                            sw2 = .true.  ! allowed 
                            exit wordloop
                        endif
                    enddo wordloop
                    if (sw2.eqv..false.) then
                        sw1 = .true. ! forbidden 
                        exit mwloop
                    endif
                enddo mwloop
                if (sw1.eqv..true.) then
                    sw0 = .true. ! forbidden 
                    exit seqloop
                endif
            enddo seqloop
            if (sw0.eqv..false.) then 
                indx(k) = 1 !.not.sw0 ! allowed
            endif
        enddo

        write (7, FMT = myprint)  (indx(k), k=1,couni)    

    enddo mainloop
    sec1 = secnds(sec)
    print *, 'The eigenvector is lies in', j-1, 'dimensional subspace; ', sec1, 'seconds to run'

    close (7)

!    deallocate(allwords)
!    deallocate(badwords)  
    deallocate(indx)
    deallocate(temp)
    deallocate(uind)
    deallocate(couind)
end subroutine minimatrix


end module dimsetup_subs         

program dimsetup 
    use dimsetup_subs
    implicit none 
    integer :: An, ex, N, total, maxmul, uni, couni, num, digs 
    integer, allocatable :: albet(:), badwds(:,:), allwds(:,:)
    integer :: k0, k1, k2, k3
    integer(long) :: mm(4)
    character(100) :: wordsfile, matfile, infile, mapfile, mulfile, comulfile
    character(25) :: myprint, dimchar
    logical :: sw0, sw1, sw2, sw3
    real *4 sec, sec1

    if(command_argument_count()<1) then
        write(*,*) 'The number of the setup file is required. Stopping.' 
        stop
    else   
        call get_command_argument(1,dimchar) 
        read(dimchar,*) num
    endif 

    write (infile, '(A5I2A4)') trim("setup"), num, trim(".day")
    print *, infile

    open(12, file = trim(infile))
    read(12, *) An
    print *, An, '# of the alphabet'   
    allocate(albet(An))
    read(12, *) (albet(k0), k0 = 1,An)
    read(12, *) N
    print *, N, 'wordlength'   
    read(12, *) ex
    print *, ex, '# forbidden words'   
    allocate(badwds(ex,N+2))
    do k0 = 1,ex
        read(12, * ) badwds(k0,1), (badwds(k0,k1), k1 =2, 1+badwds(k0,1))
    enddo
    read(12, *) digs
    print *, 'Accuracy of ', digs, ' digits'   
    close(12)   

    write (wordsfile, '(A5I2A4)') trim("words"), num, trim(".dat")
    print *, wordsfile
    write (mulfile, '(A10I2A4)') trim("multmatrix"), num, trim(".dat")
    print *, mulfile
    write (comulfile, '(A12I2A4)') trim("comultmatrix"), num, trim(".dat")
    print *, comulfile
    write (mapfile, '(A13I2A4)') trim("smallmatrices"), num, trim(".dat")
    print *, mapfile
    write (matfile, '(A10I2A4)') trim("minimatrix"), num, trim(".dat")
    print *, matfile

    total = firstsearch(badwds,ex,N,albet,An,wordsfile,mapfile)
    
    allocate(allwds(total,N))

    call wordsearch(total,N,ex,badwds,wordsfile,mulfile,comulfile,allwds,maxmul,uni,couni)
    
    call minimatrix(total,N,ex,maxmul,uni,couni,allwds,badwds,mulfile,comulfile,matfile)

    print *, 'uni:', uni, 'couni:', couni

    open(7, file = "init.dat", status = 'replace')

    write(7, FMT="(7(4xI7))") num, 128, uni, couni, maxmul, total, digs

    close(7)
    deallocate(allwds)
    deallocate(badwds)

end program dimsetup





