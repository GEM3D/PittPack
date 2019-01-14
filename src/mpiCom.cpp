#include "mpiCom.h"
#include "mpi.h"
#include <iostream>
#define METHOD 0
#include <math.h>
#define DEBUG5 1

using namespace std;

MpiCom::MpiCom( int n0, int n1, int n2, int px, int py )
{
    nx = n0;
    ny = n1;
    nz = n2;

    // since depending onthe transformation we will be permutating
    // int p0,p1;

    p0 = px;
    p1 = py;

    // thos poisson solver will be used in a flow solver
    // itertively so our start and end are the same

    topology[0] = p0;
    topology[1] = p1;
    topology[2] = 1;

    MPI_Comm_dup( MPI_COMM_WORLD, &Comm );

    MPI_Comm_rank( Comm, &myRank );
    MPI_Comm_size( Comm, &comSize );

    p0 = sqrt( comSize );
    p1 = p0;

    nChunk = p0;

    if ( p0 * p0 != comSize )
    {
        cout << " wrong comSize, go check your comSize " << p0 << endl;
        exit( 0 );
    }

    //   getChunkSize();

    // Pn.allocate(chunkSize);

    //   P=new double[fullSize];
    allocateChunks();
}

MpiCom::~MpiCom()
{
    // MPI_Comm_free(&nbrComm0);
    // MPI_Comm_free(&nbrComm1);

    // delete [] nbrs;
    //
    delete[] Nbrs[0];
    delete[] Nbrs[1];
    // delete [] nbrsXY;
}

//
// start with the zdirection and rotate and come back to it
//
void MpiCom::allocateChunks()
{
    int n[3] = {nx / p0, ny / p1, nz};

// int nChunk=p0;
#if ( DEBUG5 )
    if ( n[0] * p0 != nx )
    {
        printf( "mesh size not divisible to given chunks" );
        exit( 0 );
    }
#endif

    P.allocate( n, nChunk );

    P.getChunkSize();

#if ( METHOD == 0 )
    {
        R.allocate( n, nChunk );
#if ( DEBUG5 )
        cout << "regular send recieve-- full recv buffer is used" << endl;
#endif
    }
#else
    {
        R.allocate( n, 1 );
#if ( DEBUG5 )
        cout << "pairwise exchange-- only chunksize buffer is used" << endl;
#endif
    }
#endif

#if ( DEBUG5 )
    cout << "Total chunkSize=" << P.getChunkSize() << endl;
#endif
    P.setDirection( 2 );
}

void MpiCom::getChunkSize()
{
    int size = nx * ny * nz;

    if ( size % ( p0 * p1 ) == 0 )
    {
        size = ( nx * ny * nz ) / ( p0 * p1 );

        cout << "chunk size " << size << endl;
    }

    else
    {
        cout << "warning not divisible " << endl;
    }

    fullSize = size;

    chunkSize = size / p0;
}

void MpiCom::getNoOfChunk()
{
    // number of chunks that every processor has
    nChunk = p0;
}

void MpiCom::constructConnectivity()
{
    Nbrs[0] = new int[p0];
    Nbrs[1] = new int[p0];

    // nbrs=new int[p0];

    int id = myRank / p0;
    int id1;
    int counter = 0;

    for ( int i = 0; i < p0; i++ )
    {
        // nbrs[i]=id*p0+i;

        Nbrs[0][i] = id * p0 + i;
    }

#if ( DEBUG5 )
    for ( int i = 0; i < p0; i++ )
    {
        // cout<<"Nbr X-dir "<<id<<" myRank "<<myRank<<" nbrs "<<nbrs[i]<<endl;
        cout << "Nbr X-dir " << id << " myRank " << myRank << " nbrs " << Nbrs[0][i] << endl;
    }
#endif

    // nbrsXY=new int[p0];

    id = myRank % p0;

    for ( int i = 0; i < p0; i++ )
    {
        // nbrsXY[i]=id+i*p0;
        Nbrs[1][i] = id + i * p0;
    }

    for ( int i = 0; i < p0; i++ )
    {
        // cout<<"Nbr-Ydir "<<id<<" myRank "<<myRank<<" nbrs with stride "<<nbrsXY[i]<<endl;
        cout << "Nbr-Ydir " << id << " myRank " << myRank << " nbrs with stride " << Nbrs[1][i] << endl;
    }
}

void MpiCom::changeOwnership()
{
    MPI_Request request[p0], request1[p0];
    MPI_Status  status;

    // uses too much memory, use pairwise exhchange instead

    for ( int i = 0; i < P.size(); i++ )
    {
        P( i ) = myRank;
        R( i ) = P( i );
    }

#if ( DEBUG5 )
    cout << "sample P " << P( 0 ) << endl;
#endif
    int id = myRank / p0;
    int id1;

    int counter = 0;
    // continue from here ....

    for ( int i = 0; i < p0; i++ )
    {
        if ( myRank != Nbrs[0][i] )
        {
            MPI_Isend( &P( 0 ) + P.chunkSize * i, P.chunkSize, MPI_DOUBLE, Nbrs[0][i], myRank, Comm, &request1[counter] );
            // cout<<" myRank "<< myRank <<" destination  "<<nbrs[i]<<" index "<<P.chunkSize*i<<endl;
            counter++;
        }
    }

    counter = 0;

    for ( int i = 0; i < p0; i++ )
    {
        if ( myRank != Nbrs[0][i] )
        {
            // MPI_Irecv( &P(0)+P.chunkSize*(nbrs[i]-p0*id), P.chunkSize , MPI_DOUBLE, nbrs[i] ,nbrs[i], Comm, &request[counter] );

            MPI_Irecv( &R( 0 ) + P.chunkSize * ( Nbrs[0][i] - p0 * id ), P.chunkSize, MPI_DOUBLE, Nbrs[0][i], Nbrs[0][i], Comm,
                       &request[counter] );

            // cout<<"--------------------------------------------------------------------- "<<endl;
            // cout<<" myRank "<< myRank <<" source  "<< nbrs[i] <<" index  "<<(nbrs[i]-id*p0)*P.chunkSize  << endl;
            // cout<<"--------------------------------------------------------------------- "<<endl;
            counter++;
        }
    }

    for ( int i = 0; i < p0 - 1; i++ )
    {
        MPI_Wait( &request[i], &status );
        MPI_Wait( &request1[i], &status );
        // perform partial transposition here
    }

#if ( DEBUG5 )
    double *rt;
    P.getAddress( rt );
    cout << "address " << rt << endl;
    // addresses match
    cout << "address " << &P( 0 ) << endl;
    cout << "chunk size " << P.getChunkSize() << endl;

    if ( myRank == 2 )
    {
        for ( int i = 0; i < P.size(); i++ )
        {
            cout << myRank << " " << R( i ) << endl;
        }
    }

#endif
}

void MpiCom::changeOwnershipPairwiseExchangeZX()
{
    MPI_Request request[p0], request1[p0], request0[p0];
    MPI_Status  status;

#if ( 1 )
#if ( DEBUG5 )
    int dest[2] = {myRank + 1, myRank - 1};
    for ( int i = 0; i < 2; i++ )
    {
        cout << " myRank  " << myRank << " will send to  " << getPeriodicRank( dest[i] ) << endl;
    }
#endif
    int counter = 0;
    // MPI_Isend( &P(0), P.chunkSize, MPI_DOUBLE, getDestinationPeriodic(myRank), myRank, Comm, &request1[0] );
    int i = myRank;
    cout << " loop= " << ( p0 + 1 ) / 2 << endl;

    for ( int j = 1; j < ( p0 ) / 2 + 1; j++ )
    {
        MPI_Isend( &P( 0 ) + P.chunkSize * getPeriodicIndex( i + j ), P.chunkSize, MPI_DOUBLE, getPeriodicRank( i + j ), myRank, Comm,
                   &request1[counter] );
        MPI_Irecv( &R( 0 ), P.chunkSize, MPI_DOUBLE, getPeriodicRank( i - j ), getPeriodicRank( i - j ), Comm, &request[counter] );
        MPI_Wait( &request[0], &status );
        // MPI_Wait( &request1[0], &status );

        // write directly to the chunk that was previously send

        counter = 0;

        MPI_Isend( &P( 0 ) + P.chunkSize * getPeriodicIndex( i - j ), P.chunkSize, MPI_DOUBLE, getPeriodicRank( i - j ), myRank, Comm,
                   &request1[counter] );
        MPI_Irecv( &P( 0 ) + P.chunkSize * getPeriodicIndex( i + j ), P.chunkSize, MPI_DOUBLE, getPeriodicRank( i + j ),
                   getPeriodicRank( i + j ), Comm, &request[counter] );
        MPI_Wait( &request[0], &status );
        // MPI_Wait( &request1[0], &status );

        // transfer from container to the desired location
        // cout<<nChunk<<endl;
        /*
        i=myRank;
        cout<<"**********************************"<<endl;
        cout<<" myRank "<<myRank<<" destination "<<getPeriodicRank(i-1) <<" tag and chunk id " <<getPeriodicIndex(i-1)<<endl;
        cout<<" myRank "<<myRank<<" source "<<getPeriodicRank(i+1) <<" chunk ID " <<getPeriodicIndex(i+1)<<endl;
        cout<<"**********************************"<<endl;
        */
        for ( int k = 0; k < P.chunkSize; k++ )
        {
            P( k + P.chunkSize * getPeriodicIndex( i - j ) ) = R( k );
        }
    }
#endif

#if ( DEBUG5 )

    // cout<<"SUCCESS " <<endl;
    for ( int i = 0; i < nChunk; i++ )
    {
        cout << "ZX my rank " << myRank << " " << P( 0 + P.chunkSize * i ) << endl;
    }
#endif
    /*

    double *rt;
    P.getAddress(rt);
    cout<<"address "<<rt<<endl;
    // addresses match
    cout<<"address "<<&P(0)<<endl;
    cout<<"chunk size "<<P.getChunkSize()<<endl;
    if(myRank==2)
    {
    for(int i=0;i<P.size();i++)
    {
    //cout<<myRank<<" "<<R(i)<<endl;
    }
    }
    */
}

void MpiCom::changeOwnershipPairwiseExchangeXY()
{
    MPI_Request request[p0], request1[p0], request0[p0];
    MPI_Status  status;

    int stride = p0;
    // uses too much memory, use pairwise exhchange instead
    //
    // network congestion

#if ( DEBUG5 )
    int dest[2] = {myRank + 1, myRank - 1};
    for ( int i = 1; i < 3; i++ )
    {
        cout << " myRank  " << myRank << " will send to  " << getPeriodicRankXY( myRank / stride + i ) << endl;
    }
#endif

#if ( 1 )
    int counter = 0;
    // MPI_Isend( &P(0), P.chunkSize, MPI_DOUBLE, getDestinationPeriodic(myRank), myRank, Comm, &request1[0] );
    int i = myRank / stride;
    cout << " loop= " << ( p0 + 1 ) / 2 << endl;

    for ( int j = 1; j < ( p0 ) / 2 + 1; j++ )
    {
        MPI_Isend( &P( 0 ) + P.chunkSize * getPeriodicIndex( i + j ), P.chunkSize, MPI_DOUBLE, getPeriodicRankXY( i + j ), myRank, Comm,
                   &request1[counter] );
        MPI_Irecv( &R( 0 ), P.chunkSize, MPI_DOUBLE, getPeriodicRankXY( i - j ), getPeriodicRankXY( i - j ), Comm, &request[counter] );
        MPI_Wait( &request[0], &status );
        MPI_Wait( &request1[0], &status );

        // write directly to the chunk that was previously send

        counter = 0;

        MPI_Isend( &P( 0 ) + P.chunkSize * getPeriodicIndex( i - j ), P.chunkSize, MPI_DOUBLE, getPeriodicRankXY( i - j ), myRank, Comm,
                   &request1[counter] );
        MPI_Irecv( &P( 0 ) + P.chunkSize * getPeriodicIndex( i + j ), P.chunkSize, MPI_DOUBLE, getPeriodicRankXY( i + j ),
                   getPeriodicRankXY( i + j ), Comm, &request[counter] );
        MPI_Wait( &request[0], &status );
        MPI_Wait( &request1[0], &status );

        // transfer from container to the desired location
        // cout<<nChunk<<endl;

        for ( int k = 0; k < P.chunkSize; k++ )
        {
            P( k + P.chunkSize * getPeriodicIndex( i - j ) ) = R( k );
        }

        /*
        i=myRank/p0;
        cout<<"**********************************"<<endl;
        cout<<" myRank "<<myRank<<" destination "<<getPeriodicRankXY(i-1) <<" tag and chunk id " <<getPeriodicIndex(i-1)<<endl;
        cout<<" myRank "<<myRank<<" source "<<getPeriodicRankXY(i+1) <<" chunk ID " <<getPeriodicIndex(i+1)<<endl;
        cout<<"**********************************"<<endl;
        */
    }

#if ( DEBUG5 )
    cout << "SUCCESS " << endl;

    for ( int i = 0; i < nChunk; i++ )
    {
        cout << " XtoY my rank " << myRank << " " << P( 0 + P.chunkSize * i ) << endl;
    }
#endif

#endif
}

void MpiCom::initialize()
{
    for ( int i = 0; i < nChunk; i++ )
    {
        for ( int j = 0; j < P.chunkSize; j++ )
        {
            P( j + P.chunkSize * i ) = (myRank)*10 + i;
        }
    }
}

void MpiCom::changeOwnershipPairwiseExchangeYZ()
{
    MPI_Request request[p0], request1[p0], request0[p0];
    MPI_Status  status;

    int stride = p0;

#if ( DEBUG5 )
    int dest[2] = {myRank + 1, myRank - 1};
    for ( int i = 1; i < 3; i++ )
    {
        cout << " myRank  " << myRank << " will send to  " << getPeriodicRankXY( myRank + i ) << endl;
    }
#endif

#if ( 1 )
    int counter = 0;
    // MPI_Isend( &P(0), P.chunkSize, MPI_DOUBLE, getDestinationPeriodic(myRank), myRank, Comm, &request1[0] );
    int i = myRank;
    cout << " loop= " << ( p0 + 1 ) / 2 << endl;

    for ( int j = 1; j < p0 / 2 + 1; j++ )
    {
        MPI_Isend( &P( 0 ) + P.chunkSize * getPeriodicIndex( i + j ), P.chunkSize, MPI_DOUBLE, getPeriodicRank( i + j ), myRank, Comm,
                   &request1[counter] );
        MPI_Irecv( &R( 0 ), P.chunkSize, MPI_DOUBLE, getPeriodicRank( i - j ), getPeriodicRank( i - j ), Comm, &request[counter] );
        MPI_Wait( &request[0], &status );
        MPI_Wait( &request1[0], &status );

        // write directly to the chunk that was previously send

        counter = 0;

        MPI_Isend( &P( 0 ) + P.chunkSize * getPeriodicIndex( i - j ), P.chunkSize, MPI_DOUBLE, getPeriodicRank( i - j ), myRank, Comm,
                   &request1[counter] );
        MPI_Irecv( &P( 0 ) + P.chunkSize * getPeriodicIndex( i + j ), P.chunkSize, MPI_DOUBLE, getPeriodicRank( i + j ),
                   getPeriodicRank( i + j ), Comm, &request[counter] );
        MPI_Wait( &request[0], &status );
        MPI_Wait( &request1[0], &status );

        // transfer from container to the desired location
        // cout<<nChunk<<endl;

        for ( int k = 0; k < P.chunkSize; k++ )
        {
            P( k + P.chunkSize * getPeriodicIndex( i - j ) ) = R( k );
        }

        /*
        i=myRank/p0;
        cout<<"**********************************"<<endl;
        cout<<" myRank "<<myRank<<" destination "<<getPeriodicRankXY(i-1) <<" tag and chunk id " <<getPeriodicIndex(i-1)<<endl;
        cout<<" myRank "<<myRank<<" source "<<getPeriodicRankXY(i+1) <<" chunk ID " <<getPeriodicIndex(i+1)<<endl;
        cout<<"**********************************"<<endl;
        */
    }

#if ( 1 )
    cout << "SUCCESS " << endl;

    for ( int i = 0; i < nChunk; i++ )
    {
        cout << "my rank " << myRank << " " << P( 0 + P.chunkSize * i ) << endl;
    }
#endif

#endif
}

int MpiCom::getPeriodicRank( int rank )
{
    //
    // return the rank in a periodic way
    // 0,1,2,3,4 as the ranks, if the index is 5, 0,1,2,3,4,0,1,2,3
    //

    int index = rank;
    int size  = p0;

    // int qt=index/comSize;
    int qt = index / size;

    if ( index >= size )
    {
        index = index - qt * size;
    }
    else if ( index < 0 )
    {
        qt    = ( index + 1 ) / size;
        index = index + ( -qt + 1 ) * size;
    }

    rank = index;

    // cout<<" myrank "<<myRank<<" new index "<<rank<<" nbr[index] "<<nbrs[rank]<<endl;
    return ( Nbrs[0][rank] );
}

int MpiCom::getPeriodicRankXY( int rank )
{
    //
    // return the rank in a periodic way
    // 0,1,2,3,4 as the ranks, if the index is 5, 0,1,2,3,4,0,1,2,3
    //

    int index = rank;
    int size  = p0;

    // int qt=index/comSize;
    int qt = index / size;

    if ( index >= size )
    {
        index = index - qt * size;
    }
    else if ( index < 0 )
    {
        qt    = ( index + 1 ) / size;
        index = index + ( -qt + 1 ) * size;
    }

    rank = index;

    // cout<<" myrank "<<myRank<<" new index "<<rank<<" nbr[index] "<<nbrs[rank]<<endl;
    return ( Nbrs[1][rank] );
}

int MpiCom::getPeriodicIndex( int rank )
{
    //
    // return the rank in a periodic way
    // 0,1,2,3,4 as the ranks, if the index is 5, 0,1,2,3,4,0,1,2,3
    //

    int index = rank;
    int size  = p0;

    // int qt=index/comSize;
    int qt = index / size;

    if ( index >= size )
    {
        index = index - qt * size;
    }
    else if ( index < 0 )
    {
        qt    = ( index + 1 ) / size;
        index = index + ( -qt + 1 ) * size;
    }

    rank = index;

    // cout<<" myrank "<<myRank<<" new index "<<rank<<" nbr[index] "<<nbrs[rank]<<endl;
    return ( rank );
}

void MpiCom::graphCreate() /*!Two different communicators are required due to the presence of stride in X to Y to rotation */
{
    int indegree  = p0;
    int outdegree = p0;
    int reorder   = 0;

    // MPI_Dist_graph_create_adjacent( Comm, indegree, nbrs, MPI_UNWEIGHTED, indegree, nbrs, MPI_UNWEIGHTED, MPI_INFO_NULL, reorder,
    // &nbrComm0 );  MPI_Dist_graph_create_adjacent( Comm, indegree, nbrsXY, MPI_UNWEIGHTED, indegree, nbrsXY, MPI_UNWEIGHTED,
    // MPI_INFO_NULL, reorder, &nbrComm1 );

    MPI_Dist_graph_create_adjacent( Comm, indegree, Nbrs[0], MPI_UNWEIGHTED, indegree, Nbrs[0], MPI_UNWEIGHTED, MPI_INFO_NULL, reorder,
                                    &nbrComm[0] );
    MPI_Dist_graph_create_adjacent( Comm, indegree, Nbrs[1], MPI_UNWEIGHTED, indegree, Nbrs[1], MPI_UNWEIGHTED, MPI_INFO_NULL, reorder,
                                    &nbrComm[1] );

    // MPI_Neighbor_alltoallv( &P(0), int sendcounts[], int sdispls[], MPI_DOUBLE, recvbuf, int recvcounts[], int rdispls[], MPI_DOUBLE
    // ,nbrComm0, MPI_Request *request);

    // MPI_Neighbor_alltoallv( sendbuf, 1, MPI_INT, recvbuf, 1, MPI_INT, graphComm );

    checkGraph( 0 );
    checkGraph( 1 );
}

void MpiCom::checkGraph( int index )
{
    int nneighbors = 0;
    int indegree   = 0;
    int outdegree  = 0;

    if ( index > 1 )
    {
        printf( "unacceptable value for i, it is either zero or one" );
        exit( 0 );
    }

    int weight;

    MPI_Dist_graph_neighbors_count( nbrComm[index], &indegree, &outdegree, &weight );

    int *neighbors = new int[indegree];
    nneighbors     = indegree;

    cout << " nneighbors is = " << nneighbors << endl;
    if ( nneighbors != p0 )
    {
        cout << " nneighbors is = " << nneighbors << endl;
        exit( 0 );
    }

    int *sources       = new int[indegree];
    int *sourceweights = new int[indegree];
    int *destinations  = new int[outdegree];
    int *destweights   = new int[outdegree];

    MPI_Dist_graph_neighbors( nbrComm[index], indegree, sources, sourceweights, outdegree, destinations, destweights );

    if ( myRank == 9 )
    {
        for ( int i = 0; i < nneighbors; i++ )
        {
            cout << " myRank " << myRank << " sources " << sources[i] << endl;
            if ( Nbrs[index][i] != sources[i] )
            {
                cout << " error in neighborhood " << endl;
                exit( 0 );
            }
        }
    }
    delete[] neighbors;
    delete[] sources;
    delete[] sourceweights;
    delete[] destinations;
    delete[] destweights;
}

void MpiCom::nbrAllToAllZX() /*!Two different communicators are required due to the presence of stride in X to Y to rotation */
{
// for All to allV

#if ( DEBUG5 )
    int *sndCnts = new int[p0];
    int *disp    = new int[p0];

    for ( int i = 0; i < p0; i++ )
    {
        if ( Nbrs[0][i] != myRank )
        {
            sndCnts[i] = P.chunkSize;
        }
        else
        {
            sndCnts[i] = 0;
        }
    }

    for ( int i = 0; i < p0; i++ )
    {
        cout << " myrank " << myRank << " sendcounts " << sndCnts[i] << endl;
    }
#endif

    // MPI_Neighbor_alltoallv(&P(0), sndCnts, int sdispls[], MPI_DOUBLE, recvbuf, sndCnts, int rdispls[], MPI_DOUBLE ,nbrComm0, MPI_Request
    // *request);

    MPI_Neighbor_alltoall( &P( 0 ), P.chunkSize, MPI_DOUBLE, &R( 0 ), P.chunkSize, MPI_DOUBLE, nbrComm[0] );

    for ( int i = 0; i < nChunk * P.chunkSize; i++ )
    {
        P( i ) = R( i );
    }

#if ( 1 )
    for ( int i = 0; i < nChunk; i++ )
    {
        cout << "ZX rank " << myRank << " after transformation " << R( 0 + P.chunkSize * i ) << " " << P( 0 + P.chunkSize * i ) << endl;
    }
#endif
}

void MpiCom::nbrAllToAllXY() /*!Two different communicators are required due to the presence of stride in X to Y to rotation */
{
// for All to allV
#if ( DEBUG5 )
    int *sndCnts = new int[p0];
    int *disp    = new int[p0];

    for ( int i = 0; i < p0; i++ )
    {
        if ( Nbrs[0][i] != myRank )
        {
            sndCnts[i] = P.chunkSize;
        }
        else
        {
            sndCnts[i] = 0;
        }
    }

    for ( int i = 0; i < p0; i++ )
    {
        cout << " myrank " << myRank << " sendcounts " << sndCnts[i] << endl;
    }
#endif

    for ( int i = 0; i < nChunk; i++ )
    {
        cout << " second my rank " << myRank << " before transformation " << R( 0 + P.chunkSize * i ) << " " << P( 0 + P.chunkSize * i )
             << endl;
    }

    // MPI_Neighbor_alltoallv(&P(0), sndCnts, int sdispls[], MPI_DOUBLE, recvbuf, sndCnts, int rdispls[], MPI_DOUBLE ,nbrComm0, MPI_Request
    // *request);

    MPI_Neighbor_alltoall( &P( 0 ), P.chunkSize, MPI_DOUBLE, &R( 0 ), P.chunkSize, MPI_DOUBLE, nbrComm[1] );

    for ( int i = 0; i < P.chunkSize * nChunk; i++ )
    {
        P( i ) = R( i );
    }

#if ( 1 )
    for ( int i = 0; i < nChunk; i++ )
    {
        cout << " XY " << myRank << " after transformation " << R( 0 + P.chunkSize * i ) << " " << P( 0 + P.chunkSize * i ) << endl;
    }
#endif
}

void MpiCom::nbrAllToAllYZ() /*!Two different communicators are required due to the presence of stride in X to Y to rotation */
{
// for All to allV
#if ( DEBUG5 )
    int *sndCnts = new int[p0];
    int *disp    = new int[p0];

    for ( int i = 0; i < p0; i++ )
    {
        if ( Nbrs[0][i] != myRank )
        {
            sndCnts[i] = P.chunkSize;
        }
        else
        {
            sndCnts[i] = 0;
        }
    }

    for ( int i = 0; i < p0; i++ )
    {
        cout << " myrank " << myRank << " sendcounts " << sndCnts[i] << endl;
    }
#endif

    // MPI_Neighbor_alltoallv(&P(0), sndCnts, int sdispls[], MPI_DOUBLE, recvbuf, sndCnts, int rdispls[], MPI_DOUBLE ,nbrComm0, MPI_Request
    // *request);

    MPI_Neighbor_alltoall( &P( 0 ), P.chunkSize, MPI_DOUBLE, &R( 0 ), P.chunkSize, MPI_DOUBLE, nbrComm0 );

    for ( int i = 0; i < nChunk * P.chunkSize; i++ )
    {
        P( i ) = R( i );
    }

    cout << " Rsize " << R.size() << " R chunk size " << R.chunkSize << endl;

#if ( 1 )
    for ( int i = 0; i < nChunk; i++ )
    {
        cout << "my rank " << myRank << " after transformation " << P( 0 + P.chunkSize * i ) << endl;
    }
#endif
}

#if ( 0 )
void MpiComm::testMpiClass( MPI_Comm mpiComm )
{
    int my_rank, comsize;
    MPI_Comm_rank( MPI_COMM_WORLD, &my_rank );
    MPI_Comm_size( MPI_COMM_WORLD, &comsize );

    int nx, ny, nz;

    nx = 12;
    ny = 12;
    nz = 12;

    int p0 = 2;
    int p1 = 2;

    /*
    cout<<"enter p in x directon "<<endl;
    cin>>p0;

    cout<<"enter p in y directon "<<endl;
    cin>>p1;
    */

    MpiCom mpi( nx, ny, nz, p0, p1 );

    // test elongated in z-direction
    // since the pencil decomposition is going to be done for momentum equations, we use z-direction as the starting point
    // keeping in mind that this is going to be integrated for flow simulation

    int n[3] = {nx / p0, ny / p1, nz};

    int nChunk = 3;

    P.allocate( n, nChunk );

    P.getChunkSize();

    cout << "chunkSize=" << P.getChunkSize() << endl;

    /*
    for(int i=0;i<3;i++)
     {
    for(int j=0;j<P.getChunkSize*2;j++)
    {



     }

     }
    */
}
#endif
