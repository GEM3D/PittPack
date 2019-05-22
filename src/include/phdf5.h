#ifndef _PHDF5_H_
#define _PHDF5_H_
#include "chunkedArray.h"
#include "communicate.h"
#include "definitions.h"
#include "hdf5.h"
#include "mpi.h"
#include <string.h>

/*!
 * \class Phdf5
 *  \brief  This Writes out pencilwise decomposed data in hdf5 format in parallel with *.xmf as metadata suitable for paraview and visit
 *
 */

class Phdf5
{
    private:
    uint totalnumber;

    public:
    Phdf5(){}; /*!< Constrcutor */
    void writeMultiBlock( ChunkedArray &F, uint appx ); /*!< writes each element as block and mesh is combination of blocks, appx sets the
                                                           appendix as string for the output file  */
    void writeMultiBlockCellCenter( ChunkedArray &F, uint appx ); /*!< writes each element as block and mesh is combination of blocks, appx
                                                                     sets the appendix as
                                                                     string for the output file  */
    void xdmfMultiBlock( ChunkedArray &F, integer comsize, integer my_rank, uint offset, uint appx );
    void xdmfMultiBlockCellCenter( ChunkedArray &F, integer comsize, integer my_rank, uint offset, uint appx );
    void writeMultiBlockCellCenter( ChunkedArray &F, uint appx, int dir, int aligndir );
    void getXcoord( int *L, const double Xa, const double Xh, const int aligndir, double *xtemp );
    void getZcoord( int *L, const double Za, const double Zh, const int aligndir, double *ztemp );
    void getYcoord( int *L, const double Za, const double Zh, const int aligndir, double *ztemp );
    void getQ( ChunkedArray &F, const int chunkId, int *L, const int aligndir, double *qtemp );


   ~Phdf5(){}; /*!< Deconstructor  */
};

#endif
