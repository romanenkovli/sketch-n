!
! $Id: fpvm3.h,v 1.10 1997/07/09 13:29:35 pvmsrc Exp $
!

!  -------------------------------------------------------------------
!          PVM version 3.4:  Parallel Virtual Machine System
!                University of Tennessee, Knoxville TN.
!            Oak Ridge National Laboratory, Oak Ridge TN.
!                    Emory University, Atlanta GA.
!       Authors:  J. J. Dongarra, G. E. Fagg, M. Fischer
!           G. A. Geist, J. A. Kohl, R. J. Manchek, P. Mucci,
!          P. M. Papadopoulos, S. L. Scott, and V. S. Sunderam
!                    (C) 1997 All Rights Reserved
!
!                               NOTICE
!
!  Permission to use, copy, modify, and distribute this software and
!  its documentation for any purpose and without fee is hereby granted
!  provided that the above copyright notice appear in all copies and
!  that both the copyright notice and this permission notice appear in
!  supporting documentation.
!
!  Neither the Institutions (Emory University, Oak Ridge National
!  Laboratory, and University of Tennessee) nor the Authors make any
!  representations about the suitability of this software for any
!  purpose.  This software is provided ``as is'' without express or
!  implied warranty.
!
!  PVM version 3 was funded in part by the U.S. Department of Energy,
!  the National Science Foundation and the State of Tennessee.
!  -------------------------------------------------------------------

!     ----------------------------------
!     fpvm3.h
!
!     Definitions to be included with
!     User's Fortran application
!     ----------------------------------

      integer PVMTASKDEFAULT, PVMTASKHOST, PVMTASKARCH, PVMTASKDEBUG
      integer PVMTASKTRACE, PVMMPPFRONT, PVMHOSTCOMPL, PVMNOSPAWNPARENT
      integer PVMHOST, PVMARCH, PVMDEBUG, PVMTRACE
      integer PVMDATADEFAULT, PVMDATARAW, PVMDATAINPLACE
      integer PVMDATATRACE
      integer PVMDEFAULT, PVMRAW, PVMINPLACE
      integer PVMTASKEXIT, PVMHOSTDELETE, PVMHOSTADD, PVMROUTEADD
      integer PVMROUTEDELETE, PVMNOTIFYCANCEL
      integer PVMROUTE, PVMDEBUGMASK, PVMAUTOERR
      integer PVMOUTPUTTID, PVMOUTPUTCODE, PVMRESVTIDS
      integer PVMTRACETID, PVMTRACECODE, PVMTRACEBUFFER
      integer PVMTRACEOPTIONS, PVMFRAGSIZE, PVMSHOWTIDS, PVMNORESET
      integer PVMTRACEFULL, PVMTRACETIME, PVMTRACECOUNT
      integer PVMSOUTPUTTID, PVMSOUTPUTCODE, PVMSTRACETID
      integer PVMSTRACECODE, PVMSTRACEBUFFER, PVMSTRACEOPTIONS
      integer PVMOUTPUTCTX, PVMTRACECTX, PVMSOUTPUTCTX, PVMSTRACECTX
      integer PVMDONTROUTE, PVMALLOWDIRECT, PVMROUTEDIRECT
      integer PVMPOLLTYPE, PVMPOLLTIME, PVMPOLLCONSTANT, PVMPOLLSLEEP
      integer PVMMBOXDEFAULT, PVMMBOXPERSISTENT, PVMMBOXMULTIINSTANCE
      integer PVMMBOXOVERWRITABLE, PVMMBOXFIRSTAVAIL
      integer PVMMBOXREADANDDELETE, PVMMBOXWAITFORINFO
      integer STRING, BYTE1, INTEGER2, INTEGER4
      integer REAL4, COMPLEX8, REAL8, COMPLEX16

      integer PvmOk, PvmSysErr, PvmBadParam, PvmMismatch
      integer PvmNoData, PvmNoHost, PvmNoFile, PvmNoMem
      integer PvmBadMsg, PvmNoBuf, PvmNoSuchBuf
      integer PvmNullGroup, PvmDupGroup, PvmNoGroup
      integer PvmNotInGroup, PvmNoInst, PvmHostFail, PvmNoParent
      integer PvmNotImpl, PvmDSysErr, PvmBadVersion, PvmOutOfRes
      integer PvmDupHost, PvmCantStart, PvmAlready, PvmNoTask
      integer PvmNoEntry, PvmDupEntry, PvmOverflow, PvmDenied
      integer PvmNotFound, PvmExists, PvmHostrNMstr

!     --------------------
!     spawn 'flag' options
!     --------------------
      parameter( PVMTASKDEFAULT    =  0)
      parameter( PVMTASKHOST       =  1)
      parameter( PVMTASKARCH       =  2)
      parameter( PVMTASKDEBUG      =  4)
      parameter( PVMTASKTRACE      =  8)
      parameter( PVMMPPFRONT       = 16)
      parameter( PVMHOSTCOMPL      = 32)
      parameter( PVMNOSPAWNPARENT  = 64)

!     --------------------------------
!     old option names still supported
!     --------------------------------
      parameter( PVMHOST  =  1)
      parameter( PVMARCH  =  2)
      parameter( PVMDEBUG =  4)
      parameter( PVMTRACE =  8)

!     -------------------------
!     buffer 'encoding' options
!     -------------------------
      parameter( PVMDATADEFAULT = 0)
      parameter( PVMDATARAW     = 1)
      parameter( PVMDATAINPLACE = 2)
      parameter( PVMDATATRACE   = 4)

!     --------------------------------
!     old option names still supported
!     --------------------------------
      parameter( PVMDEFAULT = 0)
      parameter( PVMRAW     = 1)
      parameter( PVMINPLACE = 2)

!     ----------------------
!     notify 'about' options
!     ----------------------
      parameter( PVMTASKEXIT     = 1 )
      parameter( PVMHOSTDELETE   = 2 )
      parameter( PVMHOSTADD      = 3 )
      parameter( PVMROUTEADD     = 4 )
      parameter( PVMROUTEDELETE  = 5 )
      parameter( PVMNOTIFYCANCEL = 256 )

!     --------------------------------
!     packing/unpacking 'what' options
!     --------------------------------
      parameter( STRING   = 0)
      parameter( BYTE1    = 1)
      parameter( INTEGER2 = 2)
      parameter( INTEGER4 = 3)
      parameter( REAL4    = 4)
      parameter( COMPLEX8 = 5)
      parameter( REAL8    = 6)
      parameter( COMPLEX16= 7)

!     --------------------------------
!     setopt/getopt options for 'what'
!     --------------------------------
      parameter( PVMROUTE         = 1)
      parameter( PVMDEBUGMASK     = 2)
      parameter( PVMAUTOERR       = 3)
      parameter( PVMOUTPUTTID     = 4)
      parameter( PVMOUTPUTCODE    = 5)
      parameter( PVMTRACETID      = 6)
      parameter( PVMTRACECODE     = 7)
      parameter( PVMTRACEBUFFER   = 8)
      parameter( PVMTRACEOPTIONS  = 9)
      parameter( PVMFRAGSIZE      = 10)
      parameter( PVMRESVTIDS      = 11)
      parameter( PVMSOUTPUTTID    = 12)
      parameter( PVMSOUTPUTCODE   = 13)
      parameter( PVMSTRACETID     = 14)
      parameter( PVMSTRACECODE    = 15)
      parameter( PVMSTRACEBUFFER  = 16)
      parameter( PVMSTRACEOPTIONS = 17)
      parameter( PVMSHOWTIDS      = 18)
      parameter( PVMPOLLTYPE      = 19)
      parameter( PVMPOLLTIME      = 20)
      parameter( PVMOUTPUTCTX     = 21)
      parameter( PVMTRACECTX      = 22)
      parameter( PVMSOUTPUTCTX    = 23)
      parameter( PVMSTRACECTX     = 24)
      parameter( PVMNORESET       = 25)

!     --------------------------------------------
!     tracing option values for setopt function
!     --------------------------------------------
      parameter( PVMTRACEFULL     = 1)
      parameter( PVMTRACETIME     = 2)
      parameter( PVMTRACECOUNT    = 3)

!     --------------------------------------------
!     poll type options for 'how' in setopt function
!     --------------------------------------------
      parameter( PVMPOLLCONSTANT = 1)
      parameter( PVMPOLLSLEEP    = 2)

!     --------------------------------------------
!     for message mailbox operations
!     --------------------------------------------
      parameter( PVMMBOXDEFAULT        =  0)
      parameter( PVMMBOXPERSISTENT     =  1)
      parameter( PVMMBOXMULTIINSTANCE  =  2)
      parameter( PVMMBOXOVERWRITABLE   =  4)
      parameter( PVMMBOXFIRSTAVAIL     =  8)
      parameter( PVMMBOXREADANDDELETE  = 16)
      parameter( PVMMBOXWAITFORINFO    = 32)

!     --------------------------------------------
!     routing options for 'how' in setopt function
!     --------------------------------------------
      parameter( PVMDONTROUTE  = 1)
      parameter( PVMALLOWDIRECT= 2)
      parameter( PVMROUTEDIRECT= 3)

!     --------------------------
!     error 'info' return values
!     --------------------------
      parameter( PvmOk         =   0)
      parameter( PvmBadParam   =  -2)
      parameter( PvmMismatch   =  -3)
      parameter( PvmOverflow   =  -4)
      parameter( PvmNoData     =  -5)
      parameter( PvmNoHost     =  -6)
      parameter( PvmNoFile     =  -7)
      parameter( PvmDenied     =  -8)
      parameter( PvmNoMem      = -10)
      parameter( PvmBadMsg     = -12)
      parameter( PvmSysErr     = -14)
      parameter( PvmNoBuf      = -15)
      parameter( PvmNoSuchBuf  = -16)
      parameter( PvmNullGroup  = -17)
      parameter( PvmDupGroup   = -18)
      parameter( PvmNoGroup    = -19)
      parameter( PvmNotInGroup = -20)
      parameter( PvmNoInst     = -21)
      parameter( PvmHostFail   = -22)
      parameter( PvmNoParent   = -23)
      parameter( PvmNotImpl    = -24)
      parameter( PvmDSysErr    = -25)
      parameter( PvmBadVersion = -26)
      parameter( PvmOutOfRes   = -27)
      parameter( PvmDupHost    = -28)
      parameter( PvmCantStart  = -29)
      parameter( PvmAlready    = -30)
      parameter( PvmNoTask     = -31)
      parameter( PvmNotFound   = -32)
      parameter( PvmExists     = -33)
      parameter( PvmHostrNMstr = -34)

!     --------------------------
!     these are going away in the next version.
!     use the replacements
!     --------------------------
      parameter( PvmNoEntry    = -32)
      parameter( PvmDupEntry   = -33)

