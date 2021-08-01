      subroutine recola_init
#ifdef RECOLA
      use recola
      use globals_rcl, only : set_EWloop_weak_rcl,set_EWloop_QED_rcl
      implicit none
      character*200  path
      common/cmnpath/path
      double precision phmass
      common/photonmasslambda/phmass
      
      double precision q2,cllacc,cllcracc      
      double precision ame,ammu,convfac,alpha,pi
      common/parameters/ame,ammu,convfac,alpha,pi
      integer irecola,k
      common/recolamode/irecola

      double precision duv,muuv
      common/uvparams/duv,muuv

      double precision mmrcl,mercl,murcl,mdrcl,msrcl,mcrcl,mbrcl,mtrcl,
     .     mtaurcl
      common/fermassesrcl/mmrcl,mercl,mtaurcl,murcl,mdrcl,msrcl,mcrcl,
     .    mbrcl,mtrcl

      double complex rclscale
      common/rescalecouplingrecola/rclscale
      if (irecola.eq.0) return
! the line below, with .false., the amplitudes are 0 with RECOLA2.... (??)
      call set_momenta_correction_rcl(.true.)

      call set_delta_uv_rcl(duv)
      call set_mu_uv_rcl(muuv)

** from messpv_EW.F
c seems to be irrelevant  call set_delta_ir_rcl(0d0,pi**2/6d0) ! delta IR is simply NOT there in mass regularization
      
c       call set_on_shell_scheme_rcl
      call set_complex_mass_scheme_rcl
      

ccc mhmh
      call set_EWloop_QED_rcl   ! to run QED corrections only
      
      mercl = ame
      mmrcl = ammu
      mtaurcl = 1.777d0         ! from amasses in sv.F


      print*,'mtau=me or mmu in recola, for pairs w/o id particles!'
      mtaurcl = ame
c      mtaurcl = ammu
      
      murcl = 0.0052d0
      mdrcl = 0.0051d0
      mcrcl = 1.2d0
      msrcl = 0.1d0
      mtrcl = 173.2d0           ! from amasses in sv.F
      mbrcl = 4.2d0

      call set_light_fermions_rcl(1.d-8)
      
      call set_pole_mass_electron_rcl(mercl)
      call unset_light_electron_rcl

      call set_pole_mass_muon_rcl(mmrcl,0.d0)
      call unset_light_muon_rcl

      call set_pole_mass_tau_rcl(mtaurcl,0.d0) ! from amasses in sv.F
      call unset_light_tau_rcl

      call set_pole_mass_up_rcl(murcl)
      call unset_light_up_rcl

      call set_pole_mass_down_rcl(mdrcl)
      call unset_light_down_rcl

      call set_pole_mass_strange_rcl(msrcl)
      call unset_light_strange_rcl

      call set_pole_mass_charm_rcl(mcrcl,0.d0)
      call unset_light_charm_rcl

      call set_pole_mass_bottom_rcl(mbrcl,0.d0)
      call unset_light_bottom_rcl
      
      call set_pole_mass_top_rcl(mtrcl,0.d0)
      call unset_light_top_rcl
      
      call use_alpha0_scheme_rcl(alpha)
c      call use_alphaz_scheme_rcl(alpha)


      call use_mass_reg_soft_rcl(phmass)

**********      
* switching off couplings...
      call switchoff_coupling3_rcl('e-','e+','Z')
      call switchoff_coupling3_rcl('mu-','mu+','Z')
      call switchoff_coupling3_rcl('e-','e+','H')
      call switchoff_coupling3_rcl('mu-','mu+','H')
      call switchoff_coupling3_rcl('p0','p0','H')
      call switchoff_coupling3_rcl('p+','p-','H')      
      call switchoff_coupling3_rcl('p0','p0','Z')
      call switchoff_coupling3_rcl('p+','p-','Z')      
      call switchoff_coupling3_rcl('p0','p0','A')
      call switchoff_coupling3_rcl('p+','p-','A')      
      call switchoff_coupling3_rcl('W+','W-','A')
      call switchoff_coupling3_rcl('e-','e+','H')
      call switchoff_coupling3_rcl('e-','e+','p0')
      call switchoff_coupling3_rcl('mu-','mu+','p0')
      call switchoff_coupling3_rcl('e-','nu_e~','W+')
      call switchoff_coupling3_rcl('e-','nu_e~','W-')
      call switchoff_coupling3_rcl('e-','nu_e~','p+')
      call switchoff_coupling3_rcl('e-','nu_e~','p-')
      call switchoff_coupling3_rcl('e+','nu_e','W+')
      call switchoff_coupling3_rcl('e+','nu_e','W-')
      call switchoff_coupling3_rcl('e+','nu_e','p+')
      call switchoff_coupling3_rcl('e+','nu_e','p-')
      call switchoff_coupling3_rcl('mu-','nu_mu~','W+')
      call switchoff_coupling3_rcl('mu-','nu_mu~','W-')
      call switchoff_coupling3_rcl('mu-','nu_mu~','p+')
      call switchoff_coupling3_rcl('mu-','nu_mu~','p-')
      call switchoff_coupling3_rcl('mu+','nu_mu','W+')
      call switchoff_coupling3_rcl('mu+','nu_mu','W-')
      call switchoff_coupling3_rcl('mu+','nu_mu','p+')
      call switchoff_coupling3_rcl('mu+','nu_mu','p-')

c OK
c      rclscale = (1d-14,0.d0)
c      call scale_coupling3_rcl(rclscale,'tau-','tau+','A')

**********
c      call set_print_level_squared_amplitude_rcl(2)
c      print*,'PRINTING IN RECOLA-COMMENT OUT!!'
c      print*,'PRINTING IN RECOLA-COMMENT OUT!!'
c      print*,'PRINTING IN RECOLA-COMMENT OUT!!'
c      print*,'PRINTING IN RECOLA-COMMENT OUT!!'
c      print*,'PRINTING IN RECOLA-COMMENT OUT!!'      
      
      call set_output_file_rcl(path(1:len_trim(path))//'recola.out')
c      call set_output_file_rcl('*')
      call set_collier_output_dir_rcl(path(1:len_trim(path))//'cll')

** mu-e- > mu-e-gg 1)
      call define_process_rcl(11, 'mu- e- -> mu- e- A A', 'LO')
** mu+e- > mu+e-gg 2)
      call define_process_rcl(12, 'mu+ e- -> mu+ e- A A', 'LO')

** mu-e- > mu-e-g NLO 3)
      call define_process_rcl(21, 'mu- e- -> mu- e- A', 'NLO')
** mu+e- > mu+e-g NLO 4)
      call define_process_rcl(22, 'mu+ e- -> mu+ e- A', 'NLO')


** mu-e- > mu-e- NLO 5)
      call define_process_rcl(101, 'mu- e- -> mu- e-', 'NLO')


*** pairs
** mu-e- > mu-e- tau+ tau- (mtau=me or mmu) 6)
      call define_process_rcl(1001, 'mu- e- -> mu- e- tau+ tau-', 'LO')
** mu+e- > mu+e- tau+ tau- (mtau=me or mmu) 7)
      call define_process_rcl(1002, 'mu+ e- -> mu+ e- tau+ tau-', 'LO')

      
** mu-e- > mu-e- e+ e- 8)
      call define_process_rcl(1003, 'mu- e- -> mu- e- e+ e-', 'LO')
** mu+e- > mu+e- e+ e- 9)
      call define_process_rcl(1004, 'mu+ e- -> mu+ e- e+ e-', 'LO')
** mu-e- > mu-e- mu+ mu- 10)
      call define_process_rcl(1005, 'mu- e- -> mu- e- mu+ mu-', 'LO')
** mu+e- > mu+e- mu+ mu- 11)
      call define_process_rcl(1006, 'mu+ e- -> mu+ e- mu+ mu-', 'LO')

      
c      call set_draw_level_branches_rcl(1)
      
      call generate_processes_rcl


c      call get_TIs_required_accuracy_rcl(cllacc)
c      call get_TIs_critical_accuracy_rcl(cllcracc)
c      print*,cllacc,cllcracc
c      stop
      
#endif
      return
      end
***************************************************************
      subroutine exitrecola
#ifdef RECOLA
      use recola
      implicit none
      integer irecola
      common/recolamode/irecola
      if (irecola.eq.0) return
      
      call reset_recola_rcl
#endif
      end subroutine exitrecola
***************************************************************
      subroutine recola_ampl(id,p1,p2,p3,p4,k1,k2,born,oneloop)
#ifdef RECOLA
      use recola
#endif
      implicit double precision (a-h,o-z)
*      dimension p1(0:3), p2(0:3), p3(0:3), p4(0:3),
      dimension p5(0:3)
      dimension pmat(0:3,1:6),A2(2)
      double precision k1(0:3),k2(0:3)
#include "invariants.h"      
      born = 0.d0
      oneloop = 0.d0
#ifdef RECOLA
      pmat(:,1) = p1
      pmat(:,2) = p2
      pmat(:,3) = p3
      pmat(:,4) = p4
      pmat(:,5) = k1
      pmat(:,6) = k2

c      print*,'recola'
c      print*,p1
c      print*,p2
c      print*,p3
c      print*,p4
c      print*,k1
c      print*,'end recola'
      
      if (id.eq.1)  call compute_process_rcl(11,pmat(:,1:6),'LO',A2)
      if (id.eq.2)  call compute_process_rcl(12,pmat(:,1:6),'LO',A2)

      if (id.eq.3)  call compute_process_rcl(21,pmat(:,1:5),'NLO',A2)
      if (id.eq.4)  call compute_process_rcl(22,pmat(:,1:5),'NLO',A2)

      if (id.eq.5)  call compute_process_rcl(101,pmat(:,1:4),'NLO',A2)

** pairs      
      if (id.eq.6)  call compute_process_rcl(1001,pmat(:,1:6),'LO',A2)
      if (id.eq.7)  call compute_process_rcl(1002,pmat(:,1:6),'LO',A2)

      if (id.eq.8)  call compute_process_rcl(1003,pmat(:,1:6),'LO',A2)
      if (id.eq.9)  call compute_process_rcl(1004,pmat(:,1:6),'LO',A2)

      if (id.eq.10) call compute_process_rcl(1005,pmat(:,1:6),'LO',A2)
      if (id.eq.11) call compute_process_rcl(1006,pmat(:,1:6),'LO',A2)
      
      born = A2(1)
      oneloop = A2(2)
#endif
      return
      end
*******************************************
      real*16 function qtridot(qp1,qp2)
      real*16 qp1(0:3),qp2(0:3)
      qtridot = qp1(1)*qp2(1)+qp1(2)*qp2(2)+qp1(3)*qp2(3)
      return
      end
