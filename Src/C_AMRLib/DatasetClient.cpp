
//
// $Id: DatasetClient.cpp,v 1.18 2001-08-01 21:50:45 lijewski Exp $
//

#include <climits>
#include <cstring>
#include <cstdlib>
#include <cstdio>

#include <sys/types.h>
#include <sys/socket.h>
#include <sys/errno.h>
#include <netinet/in.h>
#include <netdb.h>
#include <signal.h>
#include <fcntl.h>
#include <strstream.h>
#include <unistd.h>

#include "Box.H"
#include "FArrayBox.H"
#include "MultiFab.H"
#include "TagBox.H"
#include "DatasetClient.H"

const int         MaxBufSize    = 1024;
const int         PortOffset    = 5000;
const char* const DefaultFormat = "%7.5e";

static
bool
CreateSocket (int& newsocket)
{
    int                sockfd;
    struct sockaddr_in serveraddr;
    char*              serverhost = "localhost";
    struct hostent*    serverhostp;

    int GETUID_SERVER_PORT = getuid() + PortOffset;

    if ((sockfd = socket(AF_INET, SOCK_STREAM, 0)) < 0)
    {
        perror("Bad client socket create");
        return false;
    }
    //
    // Set up the socket structures.
    //
    bzero((char*)&serveraddr, sizeof(struct sockaddr_in));

    serveraddr.sin_family = AF_INET;

    if ((serverhostp = gethostbyname(serverhost)) == (struct hostent *) 0)
    {
        cerr << "gethostbyname on " << serverhost << " failed" << endl;
        return false;
    }
    u_long sAddr(serveraddr.sin_addr.s_addr);

    bcopy(serverhostp->h_addr, (char*)&sAddr, serverhostp->h_length);

    serveraddr.sin_port = htons(GETUID_SERVER_PORT);
    //
    // Connect to the server.
    //
    if (connect(sockfd, (sockaddr *)&serveraddr, sizeof(serveraddr)) < 0)
    {
        perror ("Bad client connect");
        return false;
    }

    newsocket = sockfd;

    return true;
}

static
bool
SendString (int         sockfd,
            const char* sendstring)
{
    int  count;
    char ptrbuffer[MaxBufSize];

    if (send(sockfd, sendstring, strlen(sendstring), 0) < 0)
    {
        perror("Bad client sendstring send");
        return false;
    }
    //
    // Wait for acknowledgment.
    //
    if ((count = recv(sockfd, ptrbuffer, MaxBufSize, 0)) < 0)
    {
        perror("Bad sendstring ack.");
        return false;
    }
    ptrbuffer[count] = '\0';

    return true;
}

static
bool
SendRealArray (int        sockfd,
               Real*      data[],
               int        nvar,
               const int* lodim,
               const int* hidim)
{
    int  count;
    char buffer[MaxBufSize];
    char ptrbuffer[MaxBufSize];
    Box  dataBox = Box(IntVect(lodim), IntVect(hidim));
    //
    // Send the box.
    //
    ostrstream bufferstream(buffer, sizeof(buffer));

    bufferstream << dataBox << ends;

    if (send(sockfd, buffer, strlen(buffer), 0) < 0)
    {
        perror("Bad client box send");
        return false;
    }
    //
    // Wait for acknowledgment.
    //
    if ((count = recv(sockfd, buffer, MaxBufSize, 0)) < 0)
    {
        perror("Bad box ack.");
        return false;
    }
    buffer[count] = '\0';
    //
    // Send nComp.
    //
    sprintf(buffer, "%d", nvar);

    if (send(sockfd, buffer, strlen(buffer), 0) < 0)
    {
        perror("Bad client nComp send");
        return false;
    }
    //
    // Wait for acknowledgment.
    //
    if ((count = recv(sockfd, buffer, MaxBufSize, 0)) < 0)
    {
        perror("Bad nComp ack.");
        return false;
    }
    buffer[count] = '\0';
    //
    // Send the data.
    //
    long t_long = sizeof(Real) * dataBox.numPts();

    BL_ASSERT(t_long < INT_MAX);

    int totalDataBytes = int(t_long);

    for (int dataComponent = 0; dataComponent < nvar; dataComponent++)
    {
        int   totalBytesSent               = 0;
        int   dataBytesRemaining           = totalDataBytes;
        char* dataComponentStartingAddress = (char*)data[dataComponent];

        while (totalBytesSent < totalDataBytes)
        {
            //
            // Send a chunk of data.
            //
            char* getDataHere    = dataComponentStartingAddress+totalBytesSent;
            int   dataBufferSize = dataBytesRemaining;
            if ((count = write(sockfd, getDataHere, dataBufferSize)) < 0)
            {
                perror("Bad client data send");
                return false;
            }
            totalBytesSent     += count;
            dataBytesRemaining -= count;
        }
    }
    //
    // Send the pointer.
    //
    ostrstream ptrbufferstream(ptrbuffer, sizeof(ptrbuffer));

    ptrbufferstream << data[0] << ends;

    if (send(sockfd, ptrbuffer, strlen(ptrbuffer), 0) < 0)
    {
        perror("Bad client data ptr send");
        return false;
    }
    //
    // Wait for acknowledgment.
    //
    if ((count = recv(sockfd, ptrbuffer, MaxBufSize, 0)) < 0)
    {
        perror("Bad data ptr ack.");
        return false;
    }
    ptrbuffer[count] = '\0';

    return true;
}

bool
ArrayViewFab (FArrayBox* fab)
{
    return ArrayViewFabFormatLabel(fab, DefaultFormat, "Fab");
}

bool
ArrayView (FArrayBox* fab)
{
    return ArrayViewFab(fab);
}

bool
ArrayViewFabFormatLabel (FArrayBox*  fab,
                         const char* format,
                         const char* label)
{
    int N = fab->nComp();

    if (N < 1)
    {
        cerr << "Error in ArrayView: fab nComp < 1: fab->nComp = " << N << endl;
        return false;
    }
    if (!fab->box().ok())
    {
        cerr << "Error in ArrayView: bad fab box = "
             << fab->box() << endl;
        return false;
    }

    Real** dataArray = new Real*[N];

    for (int d = 0; d < N; d++)
        dataArray[d] = fab->dataPtr(d);  

    bool returnValue = ArrayViewRealPtrArrayNVarDims(dataArray,
                                                     N,
                                                     fab->box().smallEnd().getVect(),
                                                     fab->box().bigEnd().getVect(),
                                                     format,
                                                     label);
    delete [] dataArray;

    return returnValue;
}

bool
ArrayViewMultiFab (MultiFab* mf)
{
    return ArrayViewMultiFabFormatLabel(mf,DefaultFormat,"MultiFab");
}

bool
ArrayViewMultiFabElement (MultiFab* mf,
                          int       element)
{
    return ArrayViewMultiFabElementFormatLabel(mf,
                                               element,
                                               DefaultFormat,
                                               "MultiFab element");
}

bool
ArrayViewMultiFabElementFormatLabel (MultiFab*   mf,
                                     int         element,
                                     const char* format,
                                     const char* label)
{
    if (!mf->ok())
    {
        cerr << "Error in ArrayViewMultiFabComp: MultiFab is not ok()." << endl;
        return false;
    }
    if (element < 0 || element >= mf->length())
    {
        cerr << "Error in ArrayViewMultiFabElement:  element index is not" << endl;
        cerr << "  within range of MultiFab.length()." << endl;
        cerr << "  MultiFab.length() = " << mf->length() << endl;
        cerr << "  Requested element = " << element << endl;
        return false;
    }

    return ArrayViewFabFormatLabel(&((*mf)[element]),format,label);
}



bool
ArrayViewReal (Real*      data,
               const int* lodim,
               const int* hidim)
{
    return ArrayViewRealFormatLabel(data,
                                    lodim,
                                    hidim,
                                    DefaultFormat,
                                    "Real Array");
}

bool
ArrayViewRealFormatLabel (Real*       data,
                          const int*  lodim,
                          const int*  hidim,
                          const char* format,
                          const char* label)
{
    return ArrayViewRealNVarFormatLabel(data,1,lodim,hidim,format,label);
}

bool
ArrayViewRealNVar (Real*      data,
                   int        nvar,
                   const int* lodim,
                   const int* hidim)
{
    return ArrayViewRealNVarFormatLabel(data,
                                        nvar,
                                        lodim,
                                        hidim,
                                        DefaultFormat,
                                        "Real Array");
}

bool
ArrayViewRealNVarFormatLabel (Real*       data,
                              int         nvar,
                              const int*  lodim,
                              const int*  hidim,
                              const char* format,
                              const char* label)
{
    if (data == 0)
    {
        cerr << "Error in ArrayView: data pointer == 0" << endl;
        return false;
    }
    if (nvar < 1)
    {
        cerr << "Error in ArrayView: nComp < 1: nvar = " << nvar << endl;
        return false;
    }

    Real** dataArray = new Real*[nvar];
    long   npts      = 1; 

    for (int sd = 0; sd < BL_SPACEDIM; sd++)
        npts *= (hidim[sd] - lodim[sd] + 1);

    for (int d = 0; d < nvar; d++)
    {
        char* tempCharPtr = (char*)data;
        tempCharPtr       += d * npts * sizeof(Real);
        dataArray[d]       = (Real*)tempCharPtr;
    }

    bool returnValue = ArrayViewRealPtrArrayNVarDims(dataArray,
                                                     nvar,
                                                     lodim,
                                                     hidim,
                                                     format,
                                                     label);
    delete [] dataArray;

    return returnValue;
}

#if (BL_SPACEDIM == 2)
bool
ArrayViewRealDims (Real* data,
                   int   xlo,
                   int   ylo,
                   int   xhi,
                   int   yhi)
{
    return ArrayViewRealDimsFormatLabel(data,
                                        xlo,
                                        ylo,
                                        xhi,
                                        yhi,
                                        DefaultFormat,
                                        "Real data");
}

bool
ArrayViewRealDimsFormatLabel (Real*       data,
                              int         xlo,
                              int         ylo,
                              int         xhi,
                              int         yhi,
                              const char* format,
                              const char* label)
{
    return ArrayViewRealNVarDimsFormatLabel(data,
                                            1,
                                            xlo,
                                            ylo,
                                            xhi,
                                            yhi,
                                            format,
                                            label);
}

bool
ArrayViewRealNVarDims (Real* data,
                       int   nvar,
                       int   xlo,
                       int   ylo,
                       int   xhi,
                       int   yhi)
{
    return ArrayViewRealNVarDimsFormatLabel(data,
                                            nvar,
                                            xlo,
                                            ylo,
                                            xhi,
                                            yhi,
                                            DefaultFormat,
                                            "Real data");
}

bool
ArrayViewRealNVarDimsFormatLabel (Real*       data,
                                  int         nvar,
                                  int         xlo,
                                  int         ylo,
                                  int         xhi,
                                  int         yhi,
                                  const char* format,
                                  const char* label)
{
    int lodims[BL_SPACEDIM], hidims[BL_SPACEDIM];

    if (data == 0)
    {
        cerr << "Error in ArrayView: data pointer == 0" << endl;
        return false;
    }
    if (nvar < 1)
    {
        cerr << "Error in ArrayView: nComp < 1: nvar = " << nvar << endl;
        return false;
    }
    if (xlo > xhi)
    {
        cerr << "Error in ArrayView: xlo > xhi: " << xlo << " > " << xhi << endl;
        return false;
    }
    if (ylo > yhi)
    {
        cerr << "Error in ArrayView: ylo > yhi: " << ylo << " > " << yhi << endl;
        return false;
    }
    lodims[0] = xlo;
    lodims[1] = ylo;
    hidims[0] = xhi;
    hidims[1] = yhi;

    return ArrayViewRealNVarFormatLabel(data,nvar,lodims,hidims,format,label);
}
#else
bool
ArrayViewRealDims (Real* data,
                   int   xlo,
                   int   ylo,
                   int   zlo,
                   int   xhi,
                   int   yhi,
                   int   zhi)
{
    return ArrayViewRealDimsFormatLabel(data,
                                        xlo,
                                        ylo,
                                        zlo,
                                        xhi,
                                        yhi,
                                        zhi,
                                        DefaultFormat,
                                        "Real data");
}

bool
ArrayViewRealDimsFormatLabel (Real*       data,
                              int         xlo,
                              int         ylo,
                              int         zlo,
                              int         xhi,
                              int         yhi,
                              int         zhi,
                              const char* format,
                              const char* label)
{
    return ArrayViewRealNVarDimsFormatLabel(data,
                                            1,
                                            xlo,
                                            ylo,
                                            zlo,
                                            xhi,
                                            yhi,
                                            zhi,
                                            format,
                                            label);
}

bool
ArrayViewRealNVarDims (Real* data,
                       int   nvar,
                       int   xlo,
                       int   ylo,
                       int   zlo,
                       int   xhi,
                       int   yhi,
                       int   zhi)
{
    return ArrayViewRealNVarDimsFormatLabel(data,
                                            nvar,
                                            xlo,
                                            ylo,
                                            zlo,
                                            xhi,
                                            yhi,
                                            zhi,
                                            DefaultFormat,
                                            "Real data");
}

bool
ArrayViewRealNVarDimsFormatLabel (Real*       data,
                                  int         nvar,
                                  int         xlo,
                                  int         ylo,
                                  int         zlo,
                                  int         xhi,
                                  int         yhi,
                                  int         zhi,
                                  const char* format,
                                  const char* label)
{
    int lodims[BL_SPACEDIM], hidims[BL_SPACEDIM];

    if (data == 0)
    {
        cerr << "Error in ArrayView:  data pointer == 0" << endl;
        return false;
    }
    if (nvar < 1)
    {
        cerr << "Error in ArrayView:  nComp < 1:  nvar = " << nvar << endl;
        return false;
    }
    if (xlo > xhi)
    {
        cerr << "Error in ArrayView:  xlo > xhi:  " << xlo << " > " << xhi << endl;
        return false;
    }
    if (ylo > yhi)
    {
        cerr << "Error in ArrayView:  ylo > yhi:  " << ylo << " > " << yhi << endl;
        return false;
    }
    if (zlo > zhi)
    {
        cerr << "Error in ArrayView:  zlo > zhi:  " << zlo << " > " << zhi << endl;
        return false;
    }
    lodims[0] = xlo;
    lodims[1] = ylo;
    lodims[2] = zlo;
    hidims[0] = xhi;
    hidims[1] = yhi;
    hidims[2] = zhi;

    return ArrayViewRealNVarFormatLabel(data,nvar,lodims,hidims,format,label);
}
#endif

bool
ArrayViewRealPtrArrayNVarDims (Real*       data[],
                               int         nvar,
                               const int*  lodim,
                               const int*  hidim,
                               const char* format,
                               const char* label)
{
    int sockfd;

    if (!CreateSocket(sockfd))
        return false;
    //
    // Send data label.
    //
    if (!SendString(sockfd, label))
        return false;
    //
    // Send format.
    //
    if (!SendString(sockfd, format))
        return false;
    //
    // Send isMultiFab.
    //
    if (!SendString(sockfd, "false"))
        return false;
    //
    // Send the data.
    //
    return SendRealArray(sockfd, data, nvar, lodim, hidim);

} 

bool
ArrayViewMultiFabFormatLabel (MultiFab*   mf,
                              const char* format,
                              const char* label)
{
    int  sockfd;
    char buffer[MaxBufSize];

    if (!CreateSocket(sockfd))
        return false;
    //
    // Send data label.
    //
    if (!SendString(sockfd, label))
        return false;
    //
    // Send format.
    //
    if (!SendString(sockfd, format))
        return false;
    //
    // Send isMultiFab.
    //
    if (!SendString(sockfd, "true"))
        return false;
    //
    // Send nElements.
    //
    sprintf(buffer, "%d", mf->length());

    if (!SendString(sockfd, buffer))
        return false;
    //
    // Send the data.
    //
    for (int element = 0; element < mf->length(); element++)
    {
        //
        // Construct dataArray for this element.
        //
        FArrayBox& fab       = (*mf)[element];
        int        nvar      = fab.nComp();
        Real**     dataArray = new Real * [nvar];

        for (int d = 0; d < nvar; d++)
            dataArray[d] = fab.dataPtr(d);

        if (!SendRealArray(sockfd,
                           dataArray,
                           nvar,
                           fab.box().loVect(),
                           fab.box().hiVect()))
        {
            return false;
        }

        delete [] dataArray;
    }

    return true;
}

bool
ArrayViewTagBox (TagBox* tb)
{
    const int N = tb->nComp();

    if (N < 1)
    {
        cerr << "Error in ArrayView: fab nComp < 1: fab->nComp = " << N << endl;
        return false;
    }
    if (!tb->box().ok())
    {
        cerr << "Error in ArrayView: bad fab box = "
             << tb->box() << endl;
        return false;
    }
    //
    // Create a temp fab and put the TagBox values into it.
    //
    FArrayBox debugFab(tb->box(), N);

    for (int nv = 0; nv < N; ++nv)
    {
        Real* debugFabPtr    = debugFab.dataPtr(nv);
        char* debugTagBoxPtr = tb->dataPtr(nv);

        for (int i = 0; i < tb->box().numPts() ; ++i)
            debugFabPtr[i] = debugTagBoxPtr[i];
    }

    Real** dataArray = new Real*[N];

    for (int d = 0; d < N; d++)
        dataArray[d] = debugFab.dataPtr(d);

    bool returnValue = ArrayViewRealPtrArrayNVarDims(dataArray,
                                                     N,
                                                     debugFab.box().smallEnd().getVect(),
                                                     debugFab.box().bigEnd().getVect(),
                                                     "%3.0f",
                                                     " TagBox ");
    delete [] dataArray;

    return returnValue;
}

bool
ArrayViewTagBoxArray (TagBoxArray* tba)
{
    const int N = tba->nComp();

    if (N < 1)
    {
        cerr << "Error in ArrayView: fab nComp < 1: fab->nComp = " << N << endl;
        return false;
    }
    if (!tba->ok())
    {
        cerr << "Error in ArrayView: bad TagBoxArray." << endl;
        return false;
    }
    //
    // Create a temp MultiFab and put the TagBoxArray values into it.
    //
    MultiFab debugMultiFab(tba->boxArray(), N, tba->nGrow());

    for (int nfab = 0; nfab < tba->length(); ++nfab)
    {
        FArrayBox& debugFab    = debugMultiFab[nfab];
        TagBox   & debugTagBox = (*tba)[nfab];

        for (int nv = 0; nv < N; ++nv)
        {
            Real* debugFabPtr    = debugFab.dataPtr(nv);
            char* debugTagBoxPtr = debugTagBox.dataPtr(nv);

            for (int i = 0; i < debugTagBox.box().numPts() ; ++i)
                debugFabPtr[i] = debugTagBoxPtr[i];
        }
    }

    bool returnValue = ArrayViewMultiFabFormatLabel(&debugMultiFab,
                                                    "%3.0f",
                                                    " TagBoxArray ");
    return returnValue;
}

//
// Mumber functions of ArrayViewHelperClass -- do NOT inline these.
//

ArrayViewHelperClass::ArrayViewHelperClass () {}
ArrayViewHelperClass::~ArrayViewHelperClass () {}
