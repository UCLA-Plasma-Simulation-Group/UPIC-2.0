!-----------------------------------------------------------------------
! Null MPI library
! written by viktor k. decyk, ucla
! copyright 2017, regents of the university of california
! update: february 4, 2017
!-----------------------------------------------------------------------
      subroutine MPI_INIT(ierror)
! initialize the MPI execution environment
! ierror = error indicator
! input: none, output: ierror
      implicit none
      integer ierror
      ierror = 0
      return
      end
!-----------------------------------------------------------------------
      subroutine MPI_FINALIZE(ierror)
! terminate MPI execution environment
! ierror = error indicator
! output: ierror
      implicit none
      integer ierror
      ierror = 0
      return
      end
!-----------------------------------------------------------------------
      subroutine MPI_SEND(buf,count,datatype,dest,tag,comm,ierror)
! blocking standard mode send
! buf = initial address of send buffer
! count = number of entries to send
! datatype = datatype of each entry
! dest = rank of destination
! tag = message tag
! comm = communicator (only MPI_COMM_WORLD currently supported)
! ierror = error indicator
! input: buf, count, datatype, dest, tag, comm
! output: ierror
      implicit none
      integer buf(*)
      integer count, datatype, dest, tag, comm, ierror
      ierror = 0
      return
      end
!-----------------------------------------------------------------------
      subroutine MPI_RECV(buf,count,datatype,source,tag,comm,status,    &
     &ierror)
! blocking receive
! buf = initial address of receive buffer
! count = maximum number of entries to receive
! datatype = datatype of each entry
! source = rank of source
! tag = message tag
! comm = communicator (only MPI_COMM_WORLD currently supported)
! status = return status
! ierror = error indicator
! input: count, datatype, source, tag, comm
! output: buf, status, ierror
      implicit none
      integer buf(*), status(*)
      integer count, datatype, source, tag, comm, ierror
      ierror = 0
      return
      end
!-----------------------------------------------------------------------
      subroutine MPI_ISEND(buf,count,datatype,dest,tag,comm,request,    &
     &ierror)
! start a non-blocking send
! buf = initial address of send buffer
! count = number of entries to send
! datatype = datatype of each entry
! dest = rank of destination
! tag = message tag
! comm = communicator (only MPI_COMM_WORLD currently supported)
! request = request handle
! ierror = error indicator
! input: buf, count, datatype, dest, tag, comm
! output: request, ierror
      implicit none
      integer buf(*)
      integer count, datatype, dest, tag, comm, request, ierror
      ierror = 0
      return
      end
!-----------------------------------------------------------------------
      subroutine MPI_IRECV(buf,count,datatype,source,tag,comm,request,  &
     &ierror)
! begin a non-blocking receive
! buf = initial address of receive buffer
! count = maximum number of entries to receive
! datatype = datatype of each entry
! source = rank of source
! tag = message tag
! comm = communicator (only MPI_COMM_WORLD currently supported)
! request = request handle
! ierror = error indicator
! input: count, datatype, source, tag, comm
! output: buf, request, ierror
      implicit none
      integer buf(*)
      integer count, datatype, source, tag, comm, request, ierror
      ierror = 0
      return
      end
!-----------------------------------------------------------------------
      subroutine MPI_WAIT(request,status,ierror)
! wait for an MPI send or receive to complete
! request = request handle
! status = status object
! ierror = error indicator
! input: request
! output: request, status, ierror
      implicit none
      integer status(*)
      integer request, ierror
      ierror = 0
      return
      end
!-----------------------------------------------------------------------
      subroutine MPI_GET_COUNT(status,datatype,count,ierror)
! get the number of "top level" elements
! status = return status of receive operation
! datatype = datatype of each receive buffer entry
! count = number of received entries
! ierror = error indicator
! input: status, datatype
! output: count, ierror
      implicit none
      integer status(*)
      integer datatype, count, ierror
      ierror = 0
      return
      end
!-----------------------------------------------------------------------
      subroutine MPI_INITIALIZED(flag,ierror)
! indicate whether MPI_init has been called
! flag = true if MPI_Init has been called, false otherwise
! ierror = error indicator
! output: flag, ierror
      implicit none
      logical flag
      integer ierror
      flag = .true.
      ierror = 0
      return
      end
!-----------------------------------------------------------------------
      subroutine MPI_COMM_SIZE(comm,size,ierror)
! determine the size of the group associated with a communicator
! comm = communicator (this is ignored)
! size = number of processors in the group of comm
! ierror = error indicator
! input: comm
! output: size, ierror
      implicit none
      integer comm, size, ierror
      size = 1
      ierror = 0
      return
      end
!-----------------------------------------------------------------------
      subroutine MPI_COMM_RANK(comm,rank,ierror)
! determine the rank of the calling process in the communicator
! comm = communicator (this is ignored)
! rank = rank of the calling process in group of comm
! ierror = error indicator
! input: comm
! output: rank, ierror
      implicit none
      integer comm, rank, ierror
      rank = 0
      ierror = 0
      return
      end
!-----------------------------------------------------------------------
      subroutine MPI_BCAST(buffer,count,datatype,root,comm,ierror)
! broadcast a message from root to all processes in comm
! buffer = starting address of buffer
! count = number of entries in buffer
! datatype = datatype of buffer
! root = rank of broadcast root
! comm = communicator (only MPI_COMM_WORLD currently supported)
! ierror = error indicator
! input: buffer, count, datatype, root, comm
! output: buffer, ierror
      implicit none
      integer buffer(*)
      integer count, datatype, root, comm, ierror
      ierror = 0
      return
      end
!-----------------------------------------------------------------------
      subroutine MPI_BARRIER(comm,ierror)
! blocks each process in comm until all processes have called it.
! comm = communicator (only MPI_COMM_WORLD currently supported)
! ierror = error indicator
! input: comm
! output: ierror
      implicit none
      integer comm, ierror
      ierror = 0
      return
      end
!-----------------------------------------------------------------------
      subroutine MPI_ALLREDUCE(sendbuf,recvbuf,count,datatype,op,comm,  &
     &ierror)
! applies a reduction operation to the vector sendbuf over the set of
! processes specified by comm and places result in recvbuf on all nodes
! sendbuf = address of send buffer
! recvbuf = address of receive buffer
! count = number of elements in send buffer
! datatype = datatype of elements in send buffer
! op = reduce operation (only max, min and sum currently supported)
! comm = communicator (only MPI_COMM_WORLD currently supported)
! ierror = error indicator
! input: sendbuf, count, datatype, op, root, comm
! output: recvbuf, ierror
      implicit none
      integer sendbuf(*), recvbuf(*)
      integer count, datatype, op, comm, ierror
      ierror = 0
      return
      end
!-----------------------------------------------------------------------
      subroutine MPI_ABORT(comm,errorcode,ierror)
! force all tasks on an MPI environment to terminate
! comm = communicator (only MPI_COMM_WORLD currently supported)
! errorcode = error code to return to invoking environment
! ierror = error indicator
! input: comm, errorcode
! output: ierror
      implicit none
      integer comm, errorcode, ierror
      ierror = 0
      return
      end
!-----------------------------------------------------------------------
      function MPI_WTIME()
! return an elapsed time on the calling processor in seconds
      implicit none
      double precision MPI_WTIME
      MPI_WTIME = 0.0d0
      return
      end
