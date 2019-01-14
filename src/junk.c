#if ( 0 )
#if ( OPENACC )
#pragma acc routine
#endif
inline int PencilDcmp::newIndex( const int chunkId, const int j, const int k )
{
    return ( chunkId * nz + j * nz * ny + k ); 
}

#if ( OPENACC )
#pragma acc routine
#endif
inline void PencilDcmp::decomposeNewIndex( const int id, int &chunkId, int &j, int &k ) 
{
    j = id / ny / nz;
    chunkId = ( id % ( ny * nz ) ) / nz;
    k = ( id % ( ny * nz ) ) % nz;
}

#if ( OPENACC )
#pragma acc routine
#endif
inline void PencilDcmp::saveToTmp( const int id, double *tmp )
{
    for ( int i = 0; i < 2 * nx; i++ )
    {    
        tmp[i] = P( 2 * nx * id + i ); 
    }    
}

// given an id, will find where this chunk should be placed after the rearrangement
// and will return the chunk id that currently occupies that location
// we will save this chunk and move the current id to that location

#if ( OPENACC )
#pragma acc routine
#endif
inline int PencilDcmp::getDestinationdIndex( const int id ) 
{
    int chunkId1, j1, k1, index0;

    decomposeNewIndex( id, chunkId1, j1, k1 );

    index0 = oldIndex( chunkId1, j1, k1 );

    return ( index0 );
}


void PencilDcmp::initialize()
{

    int i1 = myRank / p0;
    int i2 = myRank % p0;

    for ( int i = 0; i < nChunk; i++ )
    {    
        for ( int j = 0; j < P.chunkSize; j++ )
        {
            // P( j + P.chunkSize *i ) = (( myRank ) * 10 + i)*(( myRank ) * 10 + i);
            //  P( j + P.chunkSize *i ) = i1*p0*p1+i2*p0+i;
            P( j + P.chunkSize *i ) = myRank * 10;
            //  cout<< " rank " <<myRank <<" i1=  " << i1 << " i2 " << i2 <<" ini tval "<<P(2*j+P.chunkSize*i) <<endl ;
        }
    }    
}

#if(0)
void PencilDcmp::changeOwnership()
{

    MPI_Request request[p0], request1[p0], request0;
    MPI_Status status;

    // uses too much memory, use pairwise exhchange instead

    for ( int i = 0; i < P.size(); i++ )
    {
        P( i ) = myRank;
        R( i ) = P( i );
    }

#if ( DEBUG )
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
            // cout<<" myRank "<< myRank <<" destination  "<<nbrs[i]<<" index
            // "<<P.chunkSize*i<<endl;
            counter++;
        }
    }

    counter = 0;

    for ( int i = 0; i < p0; i++ )
    {
        if ( myRank != Nbrs[0][i] )
        {

            // MPI_Irecv( &P(0)+P.chunkSize*(nbrs[i]-p0*id), P.chunkSize , MPI_DOUBLE,
            // nbrs[i] ,nbrs[i], Comm, &request[counter] );

            MPI_Irecv( &R( 0 ) + P.chunkSize * ( Nbrs[0][i] - p0 * id ), P.chunkSize, MPI_DOUBLE, Nbrs[0][i], Nbrs[0][i], Comm,
                       &request[counter] );

            // cout<<"---------------------------------------------------------------------
            // "<<endl;
            // cout<<" myRank "<< myRank <<" source  "<< nbrs[i] <<" index
            // "<<(nbrs[i]-id*p0)*P.chunkSize  << endl;
            // cout<<"---------------------------------------------------------------------
            // "<<endl;
            counter++;
        }
    }

 for ( int i = 0; i < p0 - 1; i++ )
    {
        MPI_Wait( &request[i], &status );
        MPI_Wait( &request1[i], &status );
        // perform partial transposition here
    }

#if ( DEBUG )
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
#endif

void PencilDcmp::printX()
{

    std::string filename = "data";
    filename.append( to_string( myRank ) );
    ofstream myfile;
    myfile.open( filename );

    // cout << "  xxxxxxx  " << myRank << " xxxxxxxxxxxxxxxxx  " << endl;
    // myfile << "  xxxxxxx  " << myRank << " xxxxxxxxxxxxxxxxx  " << endl;

    for ( int i = 0; i < nChunk * nyChunk * nzChunk; i++ )
    {
        for ( int j = 0; j < nxChunk; j++ )
        {
            //            myfile << "( " << P( 2 * ( nxChunk * i + j ) ) << " , " << P( 2 * ( nxChunk * i + j ) + 1 ) << " ) " << '\t';
            myfile << P( 2 * ( nxChunk * i + j ) ) << '\t';
        }
        // cout << endl;
        myfile << endl;
    }
    // cout << "xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx" << endl;
    myfile << "xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx" << endl;
    myfile << "xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx" << endl;
    myfile << "xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx" << endl;
    myfile << "xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx" << endl;
    myfile.close();
}
// poisson removed and made a fucntion

#if ( 0 )
#if ( COMM_PATTERN )

    M.changeOwnershipPairwiseExchangeZX();
#if ( FFTX )

    //    M.rearrange( 0, 0, 1 );
    //   M.rearrange( 0, 2 );

    MPI_Barrier( MPI_COMM_WORLD );
    t2 = MPI_Wtime();
#if ( DEBUG )
    myfile << "     Z to X rotation" << endl;
    M.printX( myfile );
#endif

// perform FFT in x direction
// remember, need to rearrange and restore each time

// step 2) change location of the array such that FFT can be performed on a contegeous array

#if ( 1 )
    M.changeLocationX();

#if ( DEBUG )
    myfile << "     change Loc in X-- where FFT should be envoked" << endl;
    M.printX( myfile );
#endif
    // step 3) perform  FFT

    M.performTransformXdir();

#if ( DEBUG )
    myfile << "     FFTX" << endl;
    M.printX( myfile );

#endif
    // step 4) restore the array to original status before FFT

    M.restoreLocationX();
#if ( DEBUG )
    myfile << "     Restore Loc in X" << endl;
    M.printX( myfile );
#endif
#endif

#if ( FFTY )
    // step 5) pencils with n(1,0,0) is converted to pencil with n(0,1,0)
    M.changeOwnershipPairwiseExchangeXY();
#if ( DEBUG )
    myfile << "     X to Y rotation" << endl;
    M.printX( myfile );
#endif
   //    M.rearrange( 0, 1 );
    // step 6) swaps X and Y coordinates, now X is Y and nx is ny

    M.rearrangeX2Y();

#if ( DEBUG )
    myfile << "     Rearrange each block from X to Y direction" << endl;

    // M.printX( myfile );
    M.printY( myfile );
#endif

    // step 7) change location of the array such that FFT can be performed on a contegeous array in the transverse direction

    M.changeLocationY();

#if ( DEBUG )
    myfile << "      change Loc in Y-- where FFT should be envoked in Y-direction" << endl;
    // M.printX( myfile );
    M.printY( myfile );

// step 8) perform  FFT transform can be performed
#endif
    M.performTransformYdir();
// M.printX( myfile );

#if ( DEBUG )
    myfile << "     FFTY" << endl;
    M.printY( myfile );

#endif

    // step 9) restore the array to original status before FFT

    M.restoreLocationY();
#if ( DEBUG )
    myfile << "      Restore Loc in Y" << endl;
    // M.printX( myfile );
    M.printY( myfile );
#endif
#endif

// step 10) pencils with n(0,1,0) is converted to pencil with n(0,0,1)
#if ( SOLVE )

    //   M.changeOwnershipPairwiseExchangeYZ();
    M.changeOwnershipPairwiseExchangeZX();

#if ( DEBUG )
    myfile << "      Rotate Y to Z " << endl;
    // M.printX( myfile );
    M.printY( myfile );
    myfile << "      EigenValue" << endl;
    M.eigenVal( myfile );
#endif
    // step 11) Customized multiBlock Thomas and periodic Thomas (with Sherman-Morrison modification)
    M.solve();
// M.printX( myfile );

#if ( DEBUG )
    myfile << "      Thomasing" << endl;
    M.printY( myfile );

#endif
    // step 12) pencils with n(0,0,1) is converted to pencil with n(0,1,0)
    // M.changeOwnershipPairwiseExchangeYZ();
    M.changeOwnershipPairwiseExchangeZX();
#if ( DEBUG )
    myfile << "      Rotate Z to Y" << endl;
    // M.printX( myfile );
    M.printY( myfile );
#endif
#endif

//****************************************************************

//                  Reverse Transform Section

//****************************************************************

#if ( IFFTY )

    // step 13) prepre for contigeuous FFT

    M.changeLocationY();
#if ( DEBUG )
    // M.printX( myfile );
    myfile << "     change Loc in Y to perform IFFTY" << endl;
    M.printY( myfile );
#endif
    // step 14) perform IFFT

    M.performInverseTransformYdir();

#if ( DEBUG )
    myfile << "     perform IFFTY" << endl;
    M.printY( myfile );
#endif

    // step 15) restore the array to original status before IFFT

    M.restoreLocationY();

#if ( DEBUG )
    myfile << "      Restore Loc in Y" << endl;
    // M.printX( myfile );
    M.printY( myfile );
#endif

    // step 16) swaps X and Y coordinates, now X is Y and nx is ny
    // need to redefine the function ??????????????????????????????????????????????????????????????????????????
    // now again switching to x-direction arrangement so printX is fine
    M.rearrangeX2YInverse();

#if ( DEBUG )
    myfile << "     Rearrange each block from Y to X direction" << endl;
    M.printX( myfile );
#endif
    // step 17) pencils with n(0,1,0) is converted to pencil with n(1,0,0)

    M.changeOwnershipPairwiseExchangeXY();
#if ( DEBUG )
    myfile << "      Rotate Y to X" << endl;
    M.printX( myfile );

#endif
#endif

#if ( IFFTX )
    // step 18) change location of the array such that IFFT can be performed on a contegeous array
    M.changeLocationX();

#if ( DEBUG )
  myfile << "     change Loc in X to perform IFFTX" << endl;

    M.printX( myfile );

#endif
    // step 19) perform  IFFT

    M.performInverseTransformXdir();

    // step 20) restore the array to original status before FFT

    M.restoreLocationX();
#if ( DEBUG )
    myfile << "     Final result" << endl;
    M.printX( myfile );
#endif
    M.rescale();

#if ( DEBUG )
    myfile << "     Final result Rescaled" << endl;
    M.printX( myfile );

#endif
#endif

#endif

#if ( OPENACC )
#pragma acc routine vector 
#endif
void PencilDcmp::changeLocationXAux(sint i)
{
#if ( OPENACC )
    double tmp[2 * NXCHUNK];
#else
    double tmp[2 * nxChunk];
#endif


        saveToTmp( jax[iax[i + 1] - 1], tmp, 0 ); 

        for ( sint j = iax[i + 1] - 1; j > iax[i]; j-- )
        {    
            // set the last one to tmp
            saveToDest( jax[j - 1], jax[j], 0 ); 
        }    
 
        saveTmpToDest( tmp, jax[iax[i]], 0 ); 
}



#if(0)
void PoissonCPU::pittPack()
{

#if ( DEBUG0 )
    ofstream myfile;

    std::string filename = "data";
    filename.append( to_string( myRank ) ); 
    //  ofstream myfile;
    myfile.open( filename );

    myfile << "xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx" << endl;
    myfile << "     Start" << endl;
    printX( myfile );

#endif

    MPI_Barrier( MPI_COMM_WORLD );
    double t1 = MPI_Wtime();

#if ( POSS )

    changeOwnershipPairwiseExchangeZX();

#if ( FFTX )
//    M.rearrange( 0, 0, 1 );
//   M.rearrange( 0, 2 );

#if ( DEBUG0 )
    myfile << "     Z to X rotation" << endl;
    printX( myfile );
#endif

    // perform FFT in x direction
    // remember, need to rearrange and restore each time

    // step 2) change location of the array such that FFT can be performed on a contegeous array

    changeLocationX();

#if ( DEBUG0 )
    myfile << "     change Loc in X-- where FFT should be envoked" << endl;
    printX( myfile );
#endif
    // step 3) perform  FFT

    performTransformXdir();

#if ( DEBUG0 )
    myfile << "     FFTX" << endl;
    printX( myfile );

#endif
    // step 4) restore the array to original status before FFT

  restoreLocationX();
#if ( DEBUG0 )
    myfile << "     Restore Loc in X" << endl;
    printX( myfile );
#endif
#endif
#if ( FFTY )
    // step 5) pencils with n(1,0,0) is converted to pencil with n(0,1,0)
    changeOwnershipPairwiseExchangeXY();
#if ( DEBUG0 )
    myfile << "     X to Y rotation" << endl;
    printX( myfile );
#endif
    //    M.rearrange( 0, 1 );
    // step 6) swaps X and Y coordinates, now X is Y and nx is ny

    rearrangeX2Y();

#if ( DEBUG0 )
    myfile << "     Rearrange each block from X to Y direction" << endl;

    // M.printX( myfile );
    printY( myfile );
#endif

    // step 7) change location of the array such that FFT can be performed on a contegeous array in the transverse direction

    changeLocationY();

#if ( DEBUG0 )
    myfile << "      change Loc in Y-- where FFT should be envoked in Y-direction" << endl;
    // M.printX( myfile );
    printY( myfile );

// step 8) perform  FFT transform can be performed
#endif
    performTransformYdir();
// M.printX( myfile );

#if ( DEBUG0 )
    myfile << "     FFTY" << endl;
    printY( myfile );

#endif

    // step 9) restore the array to original status before FFT

    restoreLocationY();
#if ( DEBUG0 )
    myfile << "      Restore Loc in Y" << endl;
    // M.printX( myfile );
    printY( myfile );
#endif

// step 10) pencils with n(0,1,0) is converted to pencil with n(0,0,1)
#if ( SOLVE )

    //   M.changeOwnershipPairwiseExchangeYZ();
    changeOwnershipPairwiseExchangeZX();

#if ( DEBUG0 )
    myfile << "      Rotate Y to Z " << endl;
    // M.printX( myfile );
    printY( myfile );
    myfile << "      EigenValue" << endl;
    eigenVal( myfile );
#endif
    // step 11) Customized multiBlock Thomas and periodic Thomas (with Sherman-Morrison modification)
    solve();
// M.printX( myfile );

#if ( DEBUG0 )
    myfile << "      Thomasing" << endl;
    printY( myfile );

#endif
    // step 12) pencils with n(0,0,1) is converted to pencil with n(0,1,0)
    // M.changeOwnershipPairwiseExchangeYZ();
    changeOwnershipPairwiseExchangeZX();
#if ( DEBUG0 )
    myfile << "      Rotate Z to Y" << endl;
    // M.printX( myfile );
    printY( myfile );
#endif
#endif

#if ( IFFTY )

    // step 13) prepre for contigeuous FFT

    changeLocationY();
#if ( DEBUG0 )
    // M.printX( myfile );
    myfile << "     change Loc in Y to perform IFFTY" << endl;
    printY( myfile );
#endif
    // step 14) perform IFFT

    performInverseTransformYdir();
  performInverseTransformYdir();

#if ( DEBUG0 )
    myfile << "     perform IFFTY" << endl;
    printY( myfile );
#endif

    // step 15) restore the array to original status before IFFT

    restoreLocationY();

#if ( DEBUG0 )
    myfile << "      Restore Loc in Y" << endl;
    // M.printX( myfile );
    printY( myfile );
#endif
    // step 16) swaps X and Y coordinates, now X is Y and nx is ny
    // need to redefine the function ??????????????????????????????????????????????????????????????????????????
    // now again switching to x-direction arrangement so printX is fine
    rearrangeX2YInverse();

#if ( DEBUG0 )
    myfile << "     Rearrange each block from Y to X direction" << endl;
    printX( myfile );
#endif
    // step 17) pencils with n(0,1,0) is converted to pencil with n(1,0,0)

    changeOwnershipPairwiseExchangeXY();
#if ( DEBUG0 )
    myfile << "      Rotate Y to X" << endl;
    printX( myfile );

#endif
#endif

#if ( IFFTX )
    // step 18) change location of the array such that IFFT can be performed on a contegeous array
    changeLocationX();

#if ( DEBUG0 )
    myfile << "     change Loc in X to perform IFFTX" << endl;

    printX( myfile );

#endif
   // step 19) perform  IFFT

    performInverseTransformXdir();

    // step 20) restore the array to original status before FFT

    restoreLocationX();
#if ( DEBUG0 )
    myfile << "     Final result" << endl;
    printX( myfile );
#endif
    rescale();

#if ( DEBUG0 )
    myfile << "     Final result Rescaled" << endl;
    printX( myfile );

#endif
#endif

    cout << " Running .......  " << endl;

#endif
#endif

#if ( DEBUG )
    myfile.close();
#endif

    MPI_Barrier( MPI_COMM_WORLD );
    double t2 = MPI_Wtime();

    if ( myRank == 0 )
    {
        cout << " Method " << COMM_PATTERN << " total_time " << t2 - t1 << endl;
    }
}

#if(0)
void PencilDcmp::solverGPU() /*!<called on CPU runs on GPU */
{

#if ( OPENACC )

    int numDevice = acc_get_num_devices( acc_device_nvidia );

    cout << endl;
    cout << " num devices " << numDevice << endl;

    if ( numDevice == 0 )
    {
        throw std::runtime_error( "Code is for Single gpu" );
    }

    cout << endl;
    cout << " num devices " << numDevice << endl;

    cout << " FFTX  " << FFTX << " FFTY " << FFTY << " IFFTX " << IFFTX << " IFFTY  " << IFFTY << endl;

#endif

#if ( DEBUG2 )
    ofstream myfile;

    std::string filename = "data";
    filename.append( to_string( myRank ) );
    //  ofstream myfile;
    myfile.open( filename );

    myfile << "xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx" << endl;
    myfile << "     Start" << endl;

    printX( myfile );

#endif

    MPI_Barrier( MPI_COMM_WORLD );
    double t1 = MPI_Wtime();

//    P.moveHostToDevice();

#pragma acc data present( P[0 : 2 * nxChunk *nyChunk *nzChunk *nChunk] )
    {

#if ( POSS )

        changeOwnershipPairwiseExchangeZX();
        P.moveHostToDevice();

#if ( FFTX )
// perform FFT in x direction
// remember, need to rearrange and restore each time

// step 2) change location of the array such that FFT can be performed on a contegeous array
#pragma acc parallel
        changeLocationX();
        // step 3) perform  FFT


// step 4) restore the array to original status before FFT
#pragma acc parallel
        restoreLocationX();

#endif

#if ( FFTY )
        // step 5) pencils with n(1,0,0) is converted to pencil with n(0,1,0)
        P.moveDeviceToHost();
        changeOwnershipPairwiseExchangeXY();
        P.moveHostToDevice();
// step 6) swaps X and Y coordinates, now X is Y and nx is ny

#pragma acc parallel
        rearrangeX2Y();

// step 7) change location of the array such that FFT can be performed on a contegeous array in the transverse direction

#pragma acc parallel
        changeLocationY();

        // step 8) perform  FFT transform can be performed
        performTransformYdir();
// M.printX( myfile );

// step 9) restore the array to original status before FFT

#pragma acc parallel
        restoreLocationY();

#endif

#if ( SOLVE )
        // step 10) pencils with n(0,1,0) is converted to pencil with n(0,0,1)

        //   M.changeOwnershipPairwiseExchangeYZ();
        P.moveDeviceToHost();
        changeOwnershipPairwiseExchangeZX();
        P.moveHostToDevice();

// step 11) Customized multiBlock Thomas and periodic Thomas (with Sherman-Morrison modification)
#pragma acc parallel
        solve();
        // M.printX( myfile );

        // step 12) pencils with n(0,0,1) is converted to pencil with n(0,1,0)
        // M.changeOwnershipPairwiseExchangeYZ();
        P.moveDeviceToHost();
        changeOwnershipPairwiseExchangeZX();
        P.moveHostToDevice();
#endif

#if ( IFFTY )
// step 13) prepre for contigeuous FFT
#pragma acc parallel
        changeLocationY();

      // step 14) perform IFFT

        performInverseTransformYdir();

// step 15) restore the array to original status before IFFT

#pragma acc parallel
        restoreLocationY();

// step 16) swaps X and Y coordinates, now X is Y and nx is ny
// need to redefine the function ??????????????????????????????????????????????????????????????????????????
// now again switching to x-direction arrangement so printX is fine
#pragma acc parallel
        rearrangeX2YInverse();

        // step 17) pencils with n(0,1,0) is converted to pencil with n(1,0,0)

        P.moveDeviceToHost();
        changeOwnershipPairwiseExchangeXY();
        P.moveHostToDevice();
#endif

#if ( IFFTX )
// step 18) change location of the array such that IFFT can be performed on a contegeous array

#pragma acc parallel
        changeLocationX();

        // step 19) perform  IFFT

        performInverseTransformXdir();

// step 20) restore the array to original status before FFT
#pragma acc parallel
        restoreLocationX();

#pragma acc parallel
        rescale();

#endif

#endif
    }

    if ( JIC )
    {
#pragma acc data present( P, R )
        {
#pragma acc parallel
            debug();
        }
    }

    P.moveDeviceToHost();

#if ( DEBUG2 )
    printX( myfile );
    myfile.close();
#endif
    cout << " *********************************" << endl;

    for ( int k = 0; k < 1; k++ )
    {
        for ( int j = 0; j < nyChunk; j++ )
        {
            for ( int i = 0; i < nxChunk; i++ )
            {
                cout << "\t " << P( i, j, k );
            }
            cout << endl;
        }
    }
}

#if ( 0 )
    dir = 0;
    M.setCoords( X, dir );
    if ( IO == 1 ) 
    {   
        M.IO( 1, dir, 0 );
    }   
// Done basically

// I do not want to rearrange the data in Z-direction, I rather customize the thomas algorithm
/*

    myfile<<"      Rearrange each block from y to z direction "<< endl;

    M.rearrange(  2, 1 );

    M.printX(myfile);


  // M.printX();
/*
    if ( my_rank == 0 )
    {
        M.printX();
    }

    //    M.changeLocationY();
    M.restoreLocationX();
    if ( my_rank == 0 )
    {
       // M.printX();
    }

    dir = 0;
    M.setCoords( X, dir );
    M.IO( 1, dir, 0 );

    dir = 1;

    M.changeOwnershipPairwiseExchangeXY();
    M.rearrange( 0, 1 );

    //    M.rearrange(1,  0, 1 );
    M.setCoords( X, dir );
    //   M.rearrange(  0, 2 );
    M.IO( 2, dir, 0 );

    if ( my_rank == 0 )
    {
        M.printY();
    }




#endif
  //    M.changeLocationY();
    M.performTransformYdir();

fstream &myfile
    M.performInverseTransformYdir();

    M.rescale();


    M.restoreLocationY();

    if ( my_rank == 0 )
    {
        M.printY();
    }

    M.IO( 3, dir, 0 );
*/
// M.getPeriodicRankXY( my_rank-p0 );
// M.getPeriodicRank( -p0/2 );
// M.getPeriodicIndexXY( my_rank+3.*p0 );

/*
    M.graphCreate();

 nt size;

    if ( direction == 0 ) 
    {   
        size = P.nx;
    }   
    else if ( direction == 1 ) 
    {   
        size = P.ny;
    }   
   int dir = 2;

    M.setCoords( X, dir );

    M.initializeTrigonometric( tags );

    M.IO( 0, dir, 0 );

    dir = 0;

    M.setCoords( X, dir );

    double t1 = MPI_Wtime();
    MPI_Barrier( MPI_COMM_WORLD );

    M.nbrAllToAllZX();

    MPI_Barrier( MPI_COMM_WORLD );
    double t2 = MPI_Wtime();

    M.IO( 1, dir, 0 );

 dir = 1;
    M.setCoords( X, dir );

    M.nbrAllToAllXY();
    M.IO( 2, dir, 0 );
*/
/*
    M.nbrAllToAllYZ();

    M.IO( 3 );
*/

/*
    if ( my_rank == 0 )
    {
        cout << " Method " << COMM_PATTERN << " total_time " << t2 - t1 << endl;
    }
*/
/*
#pragma acc data present(P)
{
#if(POSS)

//#pragma acc parallel
  M.changeOwnershipPairwiseExchangeZX();

#if ( FFTX )

//    M.rearrange( 0, 0, 1 );
//   M.rearrange( 0, 2 );

#if ( DEBUG )
    myfile << "     Z to X rotation" << endl;
    printX( myfile );
#endif

    // perform FFT in x direction
    // remember, need to rearrange and restore each time

    // step 2) change location of the array such that FFT can be performed on a contegeous array

#pragma acc parallel
    M.changeLocationX();
#endif
}
#endif
*/

PittPackResult OPENACC_Init( int &my_rank, int &com_size )
{
PittPackResult result=SUCCESS; 
#if ( OPENACC )

    acc_init( acc_device_nvidia );                                // OpenACC call

    const int num_dev = acc_get_num_devices( acc_device_nvidia ); // #GPUs

    if ( com_size == num_dev )
    {    
    const int dev_id = my_rank % num_dev;
    acc_set_device_num( dev_id, acc_device_nvidia ); // assign GPU to one MPI process
#if(DEBUG)
     cout << "MPI process " << my_rank << "  is assigned to GPU " << dev_id << "\n";
#endif
    }    
  else 
   {
      result= GPU_INIT_FAILURE;

    }    



return(result);

#endif
}



#pragma acc routine vector
void PencilDcmp::preprocessSignalDCT01( const int i, const int direction )
{

    int size;

    if ( direction == 0 )
    {
        size = nx;
    }
    else if ( direction == 1 )
    {
        size = ny;
    }

    double theta0 = pi / 2. / size;

    double theta = 0.0;

    //    Pn(0,j,k) =0.0;

    //  if ( direction == 0 )
    {
        P( 2 *i *size + 0 ) = 0.5 * cos( theta ) * P( 2 * i * size );

        P( 2 *i *size + 1 ) = 0.0;

#pragma acc loop vector

        for ( int j = 1; j < size / 2 + 1; j++ )
        {
            theta = theta0 * j;
            //
            // watch out the order, first assign the complex as I am performing in
            // place transform
            //
            /*
                       P( i, j, k, 1 ) = 0.5 * ( sin( theta ) * P( i, j, k ) - cos( theta ) * P( size - i, j, k ) );

                       P( i, j, k, 0 ) = 0.5 * ( cos( theta ) * P( i, j, k ) + sin( theta ) * P( size - i, j, k ) );
           */
            P( 2 *( i *size + j ) + 1 ) = 0.5
                                          * ( sin( theta ) * P( 2 * ( i * size + j ) ) - cos( theta ) * P( 2 * ( i * size + size - j ) ) );
            P( 2 *( i *size + j ) ) = 0.5 * ( cos( theta ) * P( 2 * ( i * size + j ) ) + sin( theta ) * P( 2 * ( i * size + size - j ) ) );
        }

#pragma acc loop vector
        for ( int j = size / 2 + 1; j < size; j++ )
        {
            //            P( i, j, k, 0 ) = P( size - i, j, k, 0 );
            //            P( i, j, k, 1 ) = -P( size - i, j, k, 1 );

            P( 2 *( i *size + j ) ) = P( 2 * ( i * size + size - j ) );
            P( 2 *( i *size + j ) + 1 ) = -P( 2 * ( i * size + size - j ) + 1 );
        }
    }
}

void PencilDcmp::postprocessSignalDCT01( const int j, const int direction )
{
    int N;

    int size;

    if ( direction == 0 )
    {
        //size = nxChunk * nChunk;
        size = nx;
    }
    else if ( direction == 1 )
    {
        //size = nyChunk * nChunk;
        size = ny;
    }

    // to prevent overwrite we use complex part and then swap real and complex
    // parts
    int upperBound = ( size - 1 ) / 2 + 1;
    int istart = upperBound;

    //   if ( direction == 0 )
    {
#pragma acc loop vector
        for ( int i = 0; i < upperBound; i++ )
        {
            // out[2 * i][0] = outC[i][0];
            //          P( 2 * i, j, k, 1 ) = P( i, j, k, 0 );
            //         P( 2 * i, j, k, 1 ) = P( i, j, k, 0 );

            P( 2 *( j *size + 2 *i ) + 1 ) = P( 2 * ( j * size + i ) );
            //            P( 2 * ( j*size+i)+ 1 ) = P( 2*(j*size+i) );
        }

        int idx = size % 2;

#pragma acc loop vector
        for ( int i = 0; i < size / 2; i++ )
        {
            //  out[2 * i + 1][0] = outC[size - i - 1][0];

            //            P( 2 * i + 1, j, k, 1 ) = P( size - i - 1, j, k, 0 );

            P( 2 *( j *size + 2 *i + 1 ) + 1 ) = P( 2 * ( j * size + size - i - 1 ) );
        }

#pragma acc loop vector
        for ( int i = 0; i < size; i++ )
        {
            // out[2 * i + 1][0] = outC[size - i - 1][0];

            // P( i, j, k, 0 ) = P( i, j, k, 1 );

            P( 2 *( i + j *size ) ) = P( 2 * ( i + j * size ) + 1 );
            // P( i, j, k, 1 ) = 0.0;

            P( 2 *( i + j *size ) + 1 ) = 0.0;
        }
  }
}


void postprocessSignalDCT10( const int i, const int direction )
{
    int size;

    if ( direction == 0 )
    {
        size = nx;
    }
    else if ( direction == 1 )
    {
        size = ny;
    }
    int idx = size % 2;

    // assign the conjugates
    if ( idx == 1 )
    {
        idx = 0;
    }
    else
    {
        idx = 1;
    }

    // make it conjugate
    //
   // if ( direction == 0 )
    {
#pragma acc loop vector
        for ( int j = 0; j < size / 2 - idx; j++ )
        {

            //           P( i + size / 2 + 1, j, k, 0 ) = P( size / 2 - i - idx, j, k, 0 );
            //           P( i + size / 2 + 1, j, k, 1 ) = -P( size / 2 - i - idx, j, k, 1 );

            P( 2 *( i *size + j + size / 2 + 1 ) ) = P( 2 * ( i * size + size / 2 - j - idx ) );
            P( 2 *( i *size + j + size / 2 + 1 ) + 1 ) = -P( 2 * ( i * size + size / 2 - j - idx ) + 1 );
        }
    }

    //   double theta0 = pi / 2 / mysize;
    double theta;

    // notice the multiplication by two here
    //    if ( direction == 0 )
    {
#pragma acc loop vector
        for ( int j = 0; j < size; j++ )
        {
            theta = ( pi / 2. / size ) * j;

            P( 2 *( i *size + j ) ) = 2. * ( P( 2 * ( i * size + j ) ) * ( cos( theta ) ) + P( 2 * ( i * size + j ) + 1 ) * sin( theta ) );
            P( 2 *( i *size + j ) + 1 ) = 0.0;
        }
    }
}
#pragma acc routine vector
void PencilDcmp::preprocessSignalDCT10( const int i, const int direction )
{

    int size;

    if ( direction == 0 )
    {
        size = nChunk * nxChunk;
    }
    else if ( direction == 1 )
    {
        size = nyChunk * nChunk;
    }

    // cout <<"szie "<<size <<endl;

    int upperBound = ( size - 1 ) / 2 + 1;
    int istart = upperBound;

#pragma acc loop vector
    for ( int j = 0; j < upperBound; j++ )
    {
        P( 2 *( i *size + j ) + 1 ) = P( 2 * ( i * size + 2 * j ) );
        //        P( 2*(i*size+ j)+1 ) = 500.;
    }

    int idx = size % 2;

#pragma acc loop vector
    for ( int j = istart; j < size; j++ )
    {
        // Pn(i,j,k) = inEx[size - 2 * ( i - istart ) - 1 - idx];
        P( 2 *( i *size + j ) + 1 ) = P( 2 * i * size + 2 * ( size - 2 * ( j - istart ) - 1 - idx ) );
        //   P(i,j,k,1)= 10000. ;
    }

    /*
        for ( int j = 0; j < size; j++ )
        {
            // Pn(i,j,k) = inEx[size - 2 * ( i - istart ) - 1 - idx];
          //  P( i, j, k, 1 ) = P(size - 2 * ( i - istart ) - 1 - idx, j, k );
         //  P(2*(i*size+j))= 10000. ;
        //   P(2*(i*size+j)+1)= 2000. ;
        }
    */

    double tmp;
#pragma acc loop vector private( tmp )
    for ( int j = 0; j < size; j++ )
  {

        tmp = P( 2 * ( i * size + j ) + 1 );
        P( 2 *( i *size + j ) ) = tmp;
        P( 2 *( i *size + j ) + 1 ) = 0.0;
    }
}




#endif


#if ( 0 )
void PoissonGPU::performTransformXdir() /*!< Called on Host and Ran on GPU*/
{
    cufftHandle plan;
    double *ptr = P.P;
#pragma acc host_data use_device( ptr )
    {

        cufftPlan1d( &plan, nx, CUFFT_Z2Z, nyChunk * nzChunk );
        cufftExecZ2Z( plan, (cufftDoubleComplex *)( ptr ), (cufftDoubleComplex *)( ptr ), CUFFT_FORWARD );
        cufftDestroy( plan );
    }
}

void PoissonGPU::performInverseTransformXdir() /*!< Called on Host and Ran on GPU*/
{
    cufftHandle plan;
    double *ptr = P.P;
#pragma acc host_data use_device( ptr )
    {

        cufftPlan1d( &plan, nx, CUFFT_Z2Z, nyChunk * nzChunk );
        cufftExecZ2Z( plan, (cufftDoubleComplex *)( ptr ), (cufftDoubleComplex *)( ptr ), CUFFT_INVERSE );
        cufftDestroy( plan );
    }
}
void PoissonGPU::performTransformYdir() /*!< Called on Host and Ran on GPU*/
{
    cufftHandle plan;
    double *ptr = P.P;
#pragma acc host_data use_device( ptr )
    {

        cufftPlan1d( &plan, ny, CUFFT_Z2Z, nxChunk * nzChunk );
        cufftExecZ2Z( plan, (cufftDoubleComplex *)( ptr ), (cufftDoubleComplex *)( ptr ), CUFFT_FORWARD );
        cufftDestroy( plan );
    }
}

void PoissonGPU::performInverseTransformYdir() /*!< Called on Host and Ran on GPU*/
{
    cufftHandle plan;
    double *ptr = P.P;
#pragma acc host_data use_device( ptr )
    {

        cufftPlan1d( &plan, ny, CUFFT_Z2Z, nxChunk * nzChunk );
        cufftExecZ2Z( plan, (cufftDoubleComplex *)( ptr ), (cufftDoubleComplex *)( ptr ), CUFFT_INVERSE );
        cufftDestroy( plan );
    }
}
#endif


#if ( 0 )
void PoissonGPU::pittPack() /*!<called on CPU runs on GPU */
{

#if ( DEBUG2 )
    ofstream myfile;

    std::string filename = "data";
    filename.append( to_string( myRank ) );
    //  ofstream myfile;
    myfile.open( filename );

    myfile << "xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx" << endl;
    myfile << "     Start" << endl;

    printX( myfile );

#endif

    MPI_Barrier( MPI_COMM_WORLD );
    double t1 = MPI_Wtime();

//    P.moveHostToDevice();

#pragma acc data present( P[0 : 2 * nxChunk *nyChunk *nzChunk *nChunk] )
    {

#if ( POSS )

        changeOwnershipPairwiseExchangeZX();
        P.moveHostToDevice();

// perform FFT in x direction
// remember, need to rearrange and restore each time

#if ( FFTX )
// step 2) change location of the array such that FFT can be performed on a contegeous array
#pragma acc parallel
        changeLocationX();
// step 3) perform  FFT

#pragma acc parallel
        preprocessSignalAccordingly( 0, 0 );

        performTransformXdir();

#pragma acc parallel
        postprocessSignalAccordingly( 0, 0 );
// step 4) restore the array to original status before FFT
#pragma acc parallel
        restoreLocationX();

#endif

#if ( FFTY )
        // step 5) pencils with n(1,0,0) is converted to pencil with n(0,1,0)
        P.moveDeviceToHost();
        changeOwnershipPairwiseExchangeXY();
        P.moveHostToDevice();
// step 6) swaps X and Y coordinates, now X is Y and nx is ny

#pragma acc parallel
        rearrangeX2Y();

// step 7) change location of the array such that FFT can be performed on a contegeous array in the transverse direction

#pragma acc parallel
        changeLocationY();

#pragma acc parallel
        preprocessSignalAccordingly( 1, 1 );

        // step 8) perform  FFT transform can be performed
        performTransformYdir();
#pragma acc parallel
        postprocessSignalAccordingly( 1, 1 );
// M.printX( myfile );

// step 9) restore the array to original status before FFT

#pragma acc parallel
        restoreLocationY();

#endif

#if ( SOLVE )
        // step 10) pencils with n(0,1,0) is converted to pencil with n(0,0,1)

        //   M.changeOwnershipPairwiseExchangeYZ();
        P.moveDeviceToHost();
        changeOwnershipPairwiseExchangeZX();
        P.moveHostToDevice();

// step 11) Customized multiBlock Thomas and periodic Thomas (with Sherman-Morrison modification)
#pragma acc parallel
        solve();
        // M.printX( myfile );

        // step 12) pencils with n(0,0,1) is converted to pencil with n(0,1,0)
        // M.changeOwnershipPairwiseExchangeYZ();
        P.moveDeviceToHost();
        changeOwnershipPairwiseExchangeZX();
        P.moveHostToDevice();
#endif

#if ( IFFTY )
// step 13) prepre for contigeuous FFT
#pragma acc parallel
        changeLocationY();

// step 14) perform IFFT

#pragma acc parallel
        preprocessSignalAccordinglyReverse( 1, 1 );

        performInverseTransformYdir();
#pragma acc parallel
        postprocessSignalAccordinglyReverse( 1, 1 );
// step 15) restore the array to original status before IFFT

#pragma acc parallel
        restoreLocationY();

// step 16) swaps X and Y coordinates, now X is Y and nx is ny
// need to redefine the function ??????????????????????????????????????????????????????????????????????????
// now again switching to x-direction arrangement so printX is fine
#pragma acc parallel
        rearrangeX2YInverse();

        // step 17) pencils with n(0,1,0) is converted to pencil with n(1,0,0)

        P.moveDeviceToHost();
        changeOwnershipPairwiseExchangeXY();
        P.moveHostToDevice();
#endif

#if ( IFFTX )
// step 18) change location of the array such that IFFT can be performed on a contegeous array

#pragma acc parallel
        changeLocationX();

// step 19) perform  IFFT

#pragma acc parallel
        preprocessSignalAccordinglyReverse( 0, 0 );

        performInverseTransformXdir();

#pragma acc parallel
        postprocessSignalAccordinglyReverse( 0, 0 );

// step 20) restore the array to original status before FFT
#pragma acc parallel
        restoreLocationX();

#pragma acc parallel
        rescale();

#endif

#endif
    }

    if ( JIC )
    {
#if ( OPENACC )
#pragma acc data present( P, R )
#endif
        {
#if ( OPENACC )
#pragma acc parallel
#endif
            debug();
        }
    }

    P.moveDeviceToHost();

#if ( DEBUG2 )
    printX( myfile );
    myfile.close();
#endif
    cout << " *********************************" << endl;
#if ( DEBUG2 )
    for ( int k = 0; k < 1; k++ )
    {
        for ( int j = 0; j < nyChunk; j++ )
        {
            for ( int i = 0; i < nxChunk; i++ )
            {
                cout << "\t " << P( i, j, k );
            }
            cout << endl;
        }
    }
#endif
}
#endif

#if ( OPENACC )
#pragma acc routine
#endif
void PencilDcmp::thomas( int i, int j, int dir, int index )
{

    double bet;

    // the enterior is always set according to eigenvalues
    // the two ends decided by BC

    double onDiag[3];
    double Sn[ZSIZE];
    /*
    #if(OPENACC)
        double Sn[ZSIZE];
    #else
        double Sn[nChunk * nzChunk];
    #endif
    */
    onDiag[1] = getEigenVal( i, j );

    //   cout << BLUE << "??????????????????????????????????my_rank " << myRank << " EigenVal " << onDiag[1] << " i= " << i + joffset
    // << " j= " << j + ioffset << endl;

    //    onDiag[0]=onDiag[1];
    //    onDiag[2]=onDiag[0];
    //    double subDiag[3]={0,1.0,1.0};
    //    double supDiag[3]{1.0,1.0,0};

    int this_rank = 1;

    //?????????????????????????????????????????????????????????????????????????????????????? need revision for Neumann
    // assign bc
    // bc 0 >> dirichlet, value known at ghost point
    if ( bc[4] == 'D' )
    {
        onDiag[0] = onDiag[1];
    }
    else if ( bc[4] == 'N' )
    {
        onDiag[0] = -1.0;
        P( 0, dir, i, j, 0, index ) = 0.0;
    }
    else if ( bc[4] == 'P' )
    {
        // peroiodic thomas (sherman-morrison) is implemented as a stand alone function
    }

    if ( bc[5] == 'D' )
    {
        onDiag[2] = onDiag[1];
    }
    else if ( bc[5] == 'N' )
    {
        onDiag[2] = -1.0;
        P( nChunk - 1, dir, i, j, nzChunk - 1, index ) = 0.0;
    }
    else if ( bc[5] == 'P' )
    {

        // peroiodic thomas (sherman-morrison) is implemented as a stand alone function
    }

#if ( DEBUG )
    if ( myRank == this_rank )
    {
        for ( int id = 0; id < nChunk; id++ )
        {
            for ( int k = 0; k < nzChunk; k++ )
            {
                cout << BLUE << "??????????????????????????????????my_rank " << myRank << " Pn [ " << 0 << "] " << P( id, dir, i, j, k, 0 )
                     << "+ " << P( id, dir, i, j, k, 1 ) << RESET << endl;
            }
        }
    }
#endif

T.thomas(P, onDiag, i, j, dir,index );

#if ( 0 )
    bet = onDiag[0];
    /*
     if(myRank==this_rank)
     {
         cout << RED << "my_rank " << myRank << " Pn [ "<< 0 <<"] " <<  P(0, i, j, 0, 0, index )  << RESET << endl;
     }
     */
    P( 0, dir, i, j, 0, index ) = P( 0, dir, i, j, 0, index ) / ( bet );

// seprate the last loop because of boundary condition

#if ( DEBUG )
    if ( myRank == this_rank )
    {
        cout << RED << "my_rank " << myRank << " Pn [ " << 0 << "] " << P( 0, dir, i, j, 0, index ) << RESET << endl;
    }
    cout << GREEN << " Nx " << nx << " Diag " << onDiag[0] << " " << onDiag[1] << " " << onDiag[2] << RESET << endl;
    cout << GREEN << " Ny " << ny << " subDiag " << subDiag[0] << " " << subDiag[1] << " " << subDiag[2] << RESET << endl;
    cout << GREEN << " Nz " << nz << " supDiag " << supDiag[0] << " " << supDiag[1] << " " << supDiag[2] << RESET << endl;
#endif
    int k = 1;
// use Sn as gamma
#if ( DEBUG )
    if ( myRank == this_rank )
    {
        cout << RED << "my_rank " << myRank << " Pn [ " << 1 << "] before " << P( 0, dir, i, j, k, index ) << RESET << endl;
    }
#endif
    Sn[k] = supDiag[0] / bet;
    bet = onDiag[1] - subDiag[1] * Sn[k];
    P( 0, dir, i, j, k, index ) = ( P( 0, dir, i, j, k, index ) - subDiag[1] * P( 0, dir, i, j, k - 1, index ) ) / bet;

#if ( DEBUG )
    if ( myRank == this_rank )
    {
        cout << RED << "my_rank " << myRank << " Pn [ " << 1 << "] " << P( 0, dir, i, j, k, index ) << RESET << endl;
    }
#endif
    /* Modify this to accommodate several chunks
        for ( int k = 2; k < n - 1; k++ )
        {
            Sn[k] = supDiag[1] / bet;
            bet = onDiag[1] - subDiag[1] * Sn[k];
            P( i, j, k, dir, index ) = ( P( i, j, k, dir, index ) - subDiag[1] * P( i, j, k - 1, dir, index ) ) / bet;
        }
    */

    int kstart = 2;
    int kend = nzChunk;
    int count = 2;

    for ( int id = 0; id < nChunk; id++ )
    {
        if ( id != 0 )
        {
            kstart = 0;
        }
        if ( id == ( nChunk - 1 ) )
        {
            kend = nzChunk - 1;
        }

        for ( int k = kstart; k < kend; k++ )
        {

            Sn[count] = supDiag[1] / bet;
            bet = onDiag[1] - subDiag[1] * Sn[count];
            P( id, dir, i, j, k, index ) = ( P( id, dir, i, j, k, index ) - subDiag[1] * P( id, dir, i, j, k - 1, index ) ) / bet;
#if ( DEBUG )
            if ( myRank == this_rank )
            {
                cout << RED << "my_rank " << myRank << " Pn [ " << count << "] " << P( id, dir, i, j, k, index ) << RESET << endl;
            }
#endif
            count++;
        }
    }

    /*
        k = N - 1;

        Sn[k] = supDiag[1] / bet;
        bet = onDiag[2] - subDiag[2] * Sn[k];
        P( i, j, k, dir, index ) = ( P( i, j, k, dir, index ) - subDiag[2] * P( i, j, k - 1, dir, index ) ) / bet;
    */
    // assigning the last element

    k = nChunk * nzChunk - 1;
    Sn[k] = supDiag[1] / bet;
    bet = onDiag[2] - subDiag[2] * Sn[k];
    P( nChunk - 1, dir, i, j, ( nzChunk - 1 ), index )
    = ( P( nChunk - 1, dir, i, j, ( nzChunk - 1 ), index ) - subDiag[2] * P( nChunk - 1, dir, i, j, nzChunk - 2, index ) ) / bet;
#if ( DEBUG )
    if ( myRank == this_rank )
    {
        cout << RED << "my_rank " << myRank << " Pn [ " << count << "] " << P( nChunk - 1, dir, i, j, nzChunk - 1, index ) << RESET << endl;
    }

    if ( myRank == 0 )
    {
        for ( int i = 0; i < nzChunk * nChunk; i++ )
        {

            cout << CYAN << "my_rank " << myRank << " Sn [ " << i << "] " << Sn[i] << RESET << endl;
        }
    }
#endif
/*
    for ( int k = ( n - 2 ); k >= 0; k-- )
    {
        P( i, j, k, dir, index ) -= Sn[k + 1] * P( i, j, k + 1, dir, index );
    }
*/

#if ( 1 )
    count = nChunk * nzChunk - 2;

    for ( int id = nChunk - 1; id >= 0; id-- )
    {

        if ( id != 0 )
        {
            kstart = 0;
        }

        if ( id == ( nChunk - 1 ) )
        {
            kend = nzChunk - 2;
        }
        else
        {

            kend = nzChunk - 1;
        }

        for ( int k = kend; k >= kstart; k-- )
        {

            P( id, dir, i, j, k, index ) -= Sn[count + 1] * P( id, dir, i, j, k + 1, index );
#if ( DEBUG )
            if ( myRank == this_rank )
            {
                cout << BLUE << "my_rank " << myRank << " Pn[   " << count << " ]" << P( id, 0, i, j, k, index ) << RESET << endl;
            }
#endif
            count--;
        }
    }

#endif

#endif
#if ( DEBUG )
    for ( int id = 0; id < nChunk; id++ )
    {
        for ( int k = 0; k < nzChunk; k++ )
        {

            // cout<<RED<<"my_rank "<<myRank<< " nz = "<<nzChunk<<" "<<P(id,0,i,j,k,0)<<RESET<<endl;
            cout << RED << "Thomas result for rank " << myRank << " " << P( id, dir, i, j, k, index ) << RESET << endl;
        }
    }
#endif
}

#if ( OPENACC )
#pragma acc routine
#endif
void PencilDcmp::thomasSingleBlock( int i, int j, int dir, int index )
{

    double bet;
    int n;

    double Sn[ZSIZE];
    // the enterior is always set according to eigenvalues
    // the two ends decided by BC
    /*
    #if(OPENACC)
        double Sn[ZSIZE];
    #else
        double Sn[nChunk * nzChunk];
    #endif
    */
    int Nx = nxChunk;
    int Ny = nyChunk;

    double onDiag[3];
    onDiag[1] = getEigenVal( i, j );

    int this_rank = 0;

    // assign bc
    // bc 0 >> dirichlet, value known at ghost point
    if ( bc[4] == 'D' )
    {

        onDiag[0] = onDiag[1];
    }
    else if ( bc[4] == 'N' )
    {
        onDiag[0] = -1.0;
        P( 0, dir, i, j, 0, index ) = 0.0;
    }
    else if ( bc[4] == 'P' )
    {

        // not implemented shemran morrision in thi way yet
    }

    if ( bc[5] == 'D' )
    {
        onDiag[2] = onDiag[1];
    }
    else if ( bc[5] == 'N' )
    {
        onDiag[2] = -1.0;
        P( 0, dir, i, j, nzChunk - 1, index ) = 0.0;
    }
    else if ( bc[5] == 'P' )
    {

        // not implemented shemran morrision in thi way yet
    }
  
     T.thomasSingleBlock(P,onDiag,i,j,dir,index);

#if ( 0 )
#if ( DEBUG )
    if ( myRank == this_rank )
    {
        for ( int id = 0; id < nChunk; id++ )
        {
            for ( int k = 0; k < nzChunk; k++ )
            {
                cout << BLUE << "??????????????????????????????????my_rank " << myRank << " Pn [ " << 0 << "] " << P( id, dir, i, j, k, 0 )
                     << "+ " << P( id, dir, i, j, k, 1 ) << RESET << endl;
            }
        }
    }
#endif

    bet = onDiag[0];
    /*
    if(myRank==this_rank)
    {
        cout << RED << "my_rank " << myRank << " Pn [ "<< 0 <<"] " <<  P(0, i, j, 0, 0, index )  << RESET << endl;
    }
    */
    P( 0, dir, i, j, 0, index ) = P( 0, dir, i, j, 0, index ) / ( bet );

// seprate the last loop because of boundary condition
#if ( DEBUG )
    if ( myRank == this_rank )
    {
        cout << RED << "my_rank " << myRank << " Pn [ " << 0 << "] " << P( 0, dir, i, j, 0, index ) << RESET << endl;
    }
    cout << GREEN << " N " << nx << " Diag " << onDiag[0] << " " << onDiag[1] << " " << onDiag[2] << RESET << endl;
    cout << GREEN << " N " << ny << " subDiag " << subDiag[0] << " " << subDiag[1] << " " << subDiag[2] << RESET << endl;
    cout << GREEN << " N " << nz << " supDiag " << supDiag[0] << " " << supDiag[1] << " " << supDiag[2] << RESET << endl;
#endif
    int k = 1;
// use Sn as gamma
#if ( DEBUG )
    if ( myRank == this_rank )
    {
        cout << RED << "my_rank " << myRank << " Pn [ " << 1 << "] before " << P( 0, dir, i, j, k, index ) << RESET << endl;
    }
#endif
    Sn[k] = supDiag[0] / bet;
    bet = onDiag[1] - subDiag[1] * Sn[k];
    P( 0, dir, i, j, k, index ) = ( P( 0, dir, i, j, k, index ) - subDiag[1] * P( 0, dir, i, j, k - 1, index ) ) / bet;

#if ( DEBUG )
    if ( myRank == this_rank )
    {
        cout << RED << "my_rank " << myRank << " Pn [ " << 1 << "] " << P( 0, dir, i, j, k, index ) << RESET << endl;
    }
#endif
    /* Mudify this to accommodate several chunks
        for ( int k = 2; k < n - 1; k++ )
        {
            Sn[k] = supDiag[1] / bet;
            bet = onDiag[1] - subDiag[1] * Sn[k];
            P( i, j, k, dir, index ) = ( P( i, j, k, dir, index ) - subDiag[1] * P( i, j, k - 1, dir, index ) ) / bet;
        }
    */

    for ( int k = 2; k < nzChunk - 1; k++ )
    {
        Sn[k] = supDiag[1] / bet;
        bet = onDiag[1] - subDiag[1] * Sn[k];
        P( 0, dir, i, j, k, index ) = ( P( 0, dir, i, j, k, index ) - subDiag[1] * P( 0, dir, i, j, k - 1, index ) ) / bet;
    }

    k = nzChunk - 1;

    Sn[k] = supDiag[1] / bet;
    bet = onDiag[2] - subDiag[2] * Sn[k];
    P( 0, dir, i, j, k, index ) = ( P( 0, dir, i, j, k, index ) - subDiag[2] * P( 0, dir, i, j, k - 1, index ) ) / bet;

    /*
        for ( int k = ( n - 2 ); k >= 0; k-- )
        {
            P( i, j, k, dir, index ) -= Sn[k + 1] * P( i, j, k + 1, dir, index );
        }
    */
    for ( int k = ( nzChunk - 2 ); k >= 0; k-- )
    {
        P( 0, dir, i, j, k, index ) -= Sn[k + 1] * P( 0, dir, i, j, k + 1, index );
    }
/*
    for ( int id = 0; id < nChunk; id++ )
    {
        for ( int k = 0; k < nzChunk; k++ )
        {

            // cout<<RED<<"my_rank "<<myRank<< " nz = "<<nzChunk<<" "<<P(id,0,i,j,k,0)<<RESET<<endl;
            cout << RED << "Thomas result for rank " << myRank << " " << P( id, 0, i, j, k, index ) << RESET << endl;
        }
    }
*/
#if ( DEBUG )
    for ( int id = 0; id < nChunk; id++ )
    {
        for ( int k = 0; k < nzChunk; k++ )
        {
            // cout<<RED<<"my_rank "<<myRank<< " nz = "<<nzChunk<<" "<<P(id,0,i,j,k,0)<<RESET<<endl;
            cout << RED << "Thomas result for rank " << myRank << " " << P( id, dir, i, j, k, index ) << RESET << endl;
        }
    }
#endif

#endif
}

#if ( OPENACC )
#pragma acc routine
#endif
void PencilDcmp::thomasPeriodic( int i, int j, int dir, int index )
{

    double bet;
    int n;

    // the enterior is always set according to eigenvalues
    // the two ends decided by BC

    double onDiag[3];

    double alpha = 1.0;
    double beta = 1.0;
    double Sn[ZSIZE];
    double Zn[ZSIZE];

    /*
    #if(OPENACC)
        double Sn[ZSIZE];
        double Zn[ZSIZE];
    #else
        double Sn[nChunk * nzChunk];
        double Zn[nChunk * nzChunk];
    #endif
    */
    // int Nx = nxChunk * nChunk;
    // int Ny = nyChunk * nChunk;

    // onDiag[1] = -2.0 + ( -2.0 + 2. * cosine( ( i + num[0] ) * pi / ( Nx + denum[0] ) ) )
    //          + ( -2.0 + 2. * cosine( ( j + num[1] ) * pi / ( Ny + denum[1] ) ) );

    onDiag[1] = getEigenVal( i, j );

  //  T.thomasPeriodic(P,onDiag, i, j, dir,index );

#if ( 1 )
    int this_rank = 0;

    // assign bc
    // bc 0 >> dirichlet, value known at ghost point

    onDiag[0] = onDiag[1];
    onDiag[2] = onDiag[1];

    double gamma = -onDiag[0];

    onDiag[0] = onDiag[0] - gamma;
    onDiag[2] = onDiag[2] - alpha * beta / gamma;

    bet = onDiag[0];
    /*
    if(myRank==this_rank)
    {
        cout << RED << "my_rank " << myRank << " Pn [ "<< 0 <<"] " <<  P(0, i, j, 0, 0, index )  << RESET << endl;
    }
    */
    P( 0, dir, i, j, 0, index ) = P( 0, dir, i, j, 0, index ) / ( bet );

// seprate the last loop because of boundary condition
#if ( DEBUG )
    if ( myRank == this_rank )
    {
        cout << RED << "my_rank " << myRank << " Pn [ " << 0 << "] " << P( 0, dir, i, j, 0, index ) << RESET << endl;
    }
    cout << GREEN << " N " << nx << " Diag " << onDiag[0] << " " << onDiag[1] << " " << onDiag[2] << RESET << endl;
    cout << GREEN << " N " << ny << " subDiag " << subDiag[0] << " " << subDiag[1] << " " << subDiag[2] << RESET << endl;
    cout << GREEN << " N " << nz << " supDiag " << supDiag[0] << " " << supDiag[1] << " " << supDiag[2] << RESET << endl;
#endif

    int k = 1;
#if ( DEBUG )
    // use Sn as gamma
    if ( myRank == this_rank )
    {
        // cout << RED << "my_rank " << myRank << " Pn [ " << 1 << "] before " << P( 0, dir, i, j, k, index ) << RESET << endl;
    }
#endif
    Sn[k] = supDiag[0] / bet;
    bet = onDiag[1] - subDiag[1] * Sn[k];
    P( 0, dir, i, j, k, index ) = ( P( 0, dir, i, j, k, index ) - subDiag[1] * P( 0, dir, i, j, k - 1, index ) ) / bet;

#if ( DEBUG )
    if ( myRank == this_rank )
    {
        // cout << RED << "my_rank " << myRank << " Pn [ " << 1 << "] " << P( 0, dir, i, j, k, index ) << RESET << endl;
    }
/* Mudify this to accommodate several chunks
    for ( int k = 2; k < n - 1; k++ )
    {
        Sn[k] = supDiag[1] / bet;
        bet = onDiag[1] - subDiag[1] * Sn[k];
        P( i, j, k, dir, index ) = ( P( i, j, k, dir, index ) - subDiag[1] * P( i, j, k - 1, dir, index ) ) / bet;
    }
*/

#endif
    int kstart = 2;
    int kend = nzChunk;
    int count = 2;

    for ( int id = 0; id < nChunk; id++ )
    {
        if ( id != 0 )
        {
            kstart = 0;
        }
        if ( id == ( nChunk - 1 ) )
        {
            kend = nzChunk - 1;
        }

        for ( int k = kstart; k < kend; k++ )
        {

            Sn[count] = supDiag[1] / bet;
            bet = onDiag[1] - subDiag[1] * Sn[count];
            P( id, dir, i, j, k, index ) = ( P( id, dir, i, j, k, index ) - subDiag[1] * P( id, dir, i, j, k - 1, index ) ) / bet;
#if ( DEBUG )
            if ( myRank == this_rank )
            {
                //          cout << RED << "my_rank " << myRank << " Pn [ " << count << "] " << P( id, dir, i, j, k, index ) << RESET <<
                // endl;
            }
#endif
            count++;
        }
    }

    /*
        k = N - 1;

        Sn[k] = supDiag[1] / bet;
        bet = onDiag[2] - subDiag[2] * Sn[k];
        P( i, j, k, dir, index ) = ( P( i, j, k, dir, index ) - subDiag[2] * P( i, j, k - 1, dir, index ) ) / bet;
    */
    // assigning the last element

    k = nChunk * nzChunk - 1;
    Sn[k] = supDiag[1] / bet;
    bet = onDiag[2] - subDiag[2] * Sn[k];
    P( nChunk - 1, dir, i, j, ( nzChunk - 1 ), index )
    = ( P( nChunk - 1, dir, i, j, ( nzChunk - 1 ), index ) - subDiag[2] * P( nChunk - 1, dir, i, j, nzChunk - 2, index ) ) / bet;

#if ( DEBUG )
    if ( myRank == this_rank )
    {
        //      cout << RED << "my_rank " << myRank << " Pn [ " << count << "] " << P( nChunk - 1, dir, i, j, nzChunk - 1, index ) << RESET
        // << endl;
    }

    if ( myRank == 0 )
    {
        for ( int i = 0; i < nzChunk * nChunk; i++ )
        {

            //        cout << CYAN << "my_rank " << myRank << " Sn [ " << i << "] " << Sn[i] << RESET << endl;
        }
    }

#endif
/*
    for ( int k = ( n - 2 ); k >= 0; k-- )
    {
        P( i, j, k, dir, index ) -= Sn[k + 1] * P( i, j, k + 1, dir, index );
    }
*/

#if ( 1 )
    count = nChunk * nzChunk - 2;

    for ( int id = nChunk - 1; id >= 0; id-- )
    {

        if ( id != 0 )
        {
            kstart = 0;
        }

        if ( id == ( nChunk - 1 ) )
        {
            kend = nzChunk - 2;
        }
        else
        {

            kend = nzChunk - 1;
        }

        for ( int k = kend; k >= kstart; k-- )
        {

            P( id, dir, i, j, k, index ) -= Sn[count + 1] * P( id, dir, i, j, k + 1, index );

            if ( myRank == this_rank )
            {
                //                cout << BLUE << "my_rank " << myRank << " Pn[   " << count << " ]" << P( id, 0, i, j, k, 0 ) << RESET <<
                // endl;
            }
            count--;
        }
    }
#endif


    count = 0;
    for ( int id = 0; id < nChunk; id++ )
    {

        for ( int k = 0; k < nzChunk; k++ )
        {
            Zn[count] = P( id, dir, i, j, k, index );
            P( id, dir, i, j, k, index ) = 0.0;
            count++;
        }
    }

    P( 0, 0, i, j, 0, index ) = gamma;
    P( nChunk - 1, 0, i, j, nzChunk - 1, index ) = alpha;

    // solve again

    bet = onDiag[0];
    /*
    if(myRank==this_rank)
    {
        cout << RED << "my_rank " << myRank << " Pn [ "<< 0 <<"] " <<  P(0, i, j, 0, 0, index )  << RESET << endl;
    }
    */
    P( 0, dir, i, j, 0, index ) = P( 0, dir, i, j, 0, index ) / ( bet );

// seprate the last loop because of boundary condition

#if ( DEBUG )
    if ( myRank == this_rank )
    {
        cout << RED << "my_rank " << myRank << " Pn [ " << 0 << "] " << P( 0, dir, i, j, 0, index ) << RESET << endl;
    }
    cout << GREEN << " N " << nx << " Diag " << onDiag[0] << " " << onDiag[1] << " " << onDiag[2] << RESET << endl;
    cout << GREEN << " N " << ny << " subDiag " << subDiag[0] << " " << subDiag[1] << " " << subDiag[2] << RESET << endl;
    cout << GREEN << " N " << nz << " supDiag " << supDiag[0] << " " << supDiag[1] << " " << supDiag[2] << RESET << endl;

#endif
    k = 1;
// use Sn as gamma
#if ( DEBUG )
    if ( myRank == this_rank )
    {
        cout << RED << "my_rank " << myRank << " Pn [ " << 1 << "] before " << P( 0, dir, i, j, k, index ) << RESET << endl;
    }
#endif
    Sn[k] = supDiag[0] / bet;
    bet = onDiag[1] - subDiag[1] * Sn[k];
    P( 0, dir, i, j, k, index ) = ( P( 0, dir, i, j, k, index ) - subDiag[1] * P( 0, dir, i, j, k - 1, index ) ) / bet;

#if ( DEBUG )
    if ( myRank == this_rank )
    {
        cout << RED << "my_rank " << myRank << " Pn [ " << 1 << "] " << P( 0, dir, i, j, k, index ) << RESET << endl;
    }
/* Mudify this to accommodate several chunks
    for ( int k = 2; k < n - 1; k++ )
    {
        Sn[k] = supDiag[1] / bet;
        bet = onDiag[1] - subDiag[1] * Sn[k];
        P( i, j, k, dir, index ) = ( P( i, j, k, dir, index ) - subDiag[1] * P( i, j, k - 1, dir, index ) ) / bet;
    }
*/

#endif
    kstart = 2;
    kend = nzChunk;
    count = 2;

    for ( int id = 0; id < nChunk; id++ )
    {
        if ( id != 0 )
        {
            kstart = 0;
        }
        if ( id == ( nChunk - 1 ) )
        {
            kend = nzChunk - 1;
        }

        for ( int k = kstart; k < kend; k++ )
        {

            Sn[count] = supDiag[1] / bet;
            bet = onDiag[1] - subDiag[1] * Sn[count];
            P( id, dir, i, j, k, index ) = ( P( id, dir, i, j, k, index ) - subDiag[1] * P( id, dir, i, j, k - 1, index ) ) / bet;

#if ( DEBUG )
            if ( myRank == this_rank )
            {
                cout << RED << "my_rank " << myRank << " Pn [ " << count << "] " << P( id, dir, i, j, k, index ) << RESET << endl;
            }
#endif
            count++;
        }
    }

    /*
        k = N - 1;

        Sn[k] = supDiag[1] / bet;
        bet = onDiag[2] - subDiag[2] * Sn[k];
        P( i, j, k, dir, index ) = ( P( i, j, k, dir, index ) - subDiag[2] * P( i, j, k - 1, dir, index ) ) / bet;
    */
    // assigning the last element

    k = nChunk * nzChunk - 1;
    Sn[k] = supDiag[1] / bet;
    bet = onDiag[2] - subDiag[2] * Sn[k];
    P( nChunk - 1, dir, i, j, ( nzChunk - 1 ), index )
    = ( P( nChunk - 1, dir, i, j, ( nzChunk - 1 ), index ) - subDiag[2] * P( nChunk - 1, dir, i, j, nzChunk - 2, index ) ) / bet;

#if ( DEBUG )
    if ( myRank == this_rank )
    {
        cout << RED << "my_rank " << myRank << " Pn [ " << count << "] " << P( nChunk - 1, dir, i, j, nzChunk - 1, index ) << RESET << endl;
    }

    if ( myRank == 0 )
    {
        for ( int i = 0; i < nzChunk * nChunk; i++ )
        {

            cout << CYAN << "my_rank " << myRank << " Sn [ " << i << "] " << Sn[i] << RESET << endl;
        }
    }
#endif
/*
    for ( int k = ( n - 2 ); k >= 0; k-- )
    {
        P( i, j, k, dir, index ) -= Sn[k + 1] * P( i, j, k + 1, dir, index );
    }
*/

#if ( 1 )
    count = nChunk * nzChunk - 2;

    for ( int id = nChunk - 1; id >= 0; id-- )
    {

        if ( id != 0 )
        {
            kstart = 0;
        }

        if ( id == ( nChunk - 1 ) )
        {
            kend = nzChunk - 2;
        }
        else
        {

            kend = nzChunk - 1;
        }

        for ( int k = kend; k >= kstart; k-- )
        {

            P( id, dir, i, j, k, index ) -= Sn[count + 1] * P( id, dir, i, j, k + 1, index );

#if ( DEBUG )
            if ( myRank == this_rank )
            {
                cout << BLUE << "my_rank " << myRank << " Pn[   " << count << " ]" << P( id, 0, i, j, k, 0 ) << RESET << endl;
            }
#endif
            count--;
        }
    }
#endif

    double fact = ( Zn[0] + beta * Zn[nChunk * nzChunk - 1] / gamma )
                  / ( 1. + P( 0, dir, i, j, 0, index ) + beta * P( nChunk - 1, dir, i, j, nzChunk - 1, index ) / gamma );

#if ( DEBUG )
    cout << " fact " << fact << endl;
#endif
    count = 0;
    for ( int id = 0; id < nChunk; id++ )
    {
        for ( int k = 0; k < nzChunk; k++ )
        {
            P( id, dir, i, j, k, index ) = Zn[count] - fact * P( id, dir, i, j, k, index );
            count++;
        }
    }

#if ( DEBUG )

    for ( int id = 0; id < nChunk; id++ )
    {
        for ( int k = 0; k < nzChunk; k++ )
        {

            // cout<<RED<<"my_rank "<<myRank<< " nz = "<<nzChunk<<" "<<P(id,0,i,j,k,0)<<RESET<<endl;
            cout << RED << "my_rank " << myRank << " " << P( id, dir, i, j, k, index ) << RESET << endl;
        }
    }
#endif

#endif

}


void PencilDcmp::changeOwnershipPairwiseExchangeYZ()
{

    MPI_Request request[p0], request1[p0], request0[p0];
    MPI_Status status;

    int stride = p0;

#if ( DEBUG )
    int dest[2] = {myRank + 1, myRank - 1};
    for ( int i = 1; i < 3; i++ )
    {
        //       cout << " myRank  " << myRank << " will send to  " << getPeriodicRank( myRank + i,1 ) << endl;
    }
#endif

#if ( 1 )
    int counter = 0;
    // MPI_Isend( &P(0), P.chunkSize, MPI_DOUBLE, getDestinationPeriodic(myRank),
    // myRank, Comm, &request1[0] );
    int i = myRank;
    //  cout << " loop= " << ( p0 + 1 ) / 2 << endl;

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
            P( k + P.chunkSize *getPeriodicIndex( i - j ) ) = R( k );
        }

        /*
        i=myRank/p0;
        cout<<"**********************************"<<endl;
        cout<<" myRank "<<myRank<<" destination "<<getPeriodicRankXY(i-1) <<" tag
        and chunk id " <<getPeriodicIndex(i-1)<<endl;
        cout<<" myRank "<<myRank<<" source "<<getPeriodicRankXY(i+1) <<" chunk ID "
        <<getPeriodicIndex(i+1)<<endl;
        cout<<"**********************************"<<endl;
        */
    }

#endif
}


#if ( 0 )
void PencilDcmp::rearrange( const int chunkId, const int original, const int next )
{

#if ( DEBUG )
    if ( chunkId >= nChunk )
    {
        cout << " invalid chunkId in   " << endl;
        exit( 0 );
    }
#endif

    /*
          for ( int i = 0; i <  nxChunk * nyChunk * nzChunk; i++ )
          {
              tmp[2*i] = P( chunkId * chunkSize + 2*i );
              tmp[2*i+1] = P( chunkId * chunkSize + 2*i+1 );
          }
  */

    for ( int k = 0; k < nzChunk; k++ )
    {
#if ( OPENACC )
#pragma acc loop vector collapse( 2 )
#endif

        cout << "=============================" << endl;
        for ( int j = 0; j < nyChunk; j++ )
        {
            for ( int i = 0; i < nxChunk; i++ )
            {
                tmp[2 * ( nxChunk * nyChunk * k + nxChunk * j + i )] = P( chunkId, original, i, j, k, 0 );

                cout << " " << tmp[2 * ( nxChunk * nyChunk * k + nxChunk * j + i )] << '\t';
            }
            cout << endl;
        }
    }

    /*
        if ( myRank == 0 )
        {
            cout << "inside here" << endl;

            for ( int k = 0; k < nzChunk; k++ )
            {
    #if(OPENACC)
    #pragma acc loop vector collapse( 2 )
    #endif
                cout << "=============================" << endl;

                for ( int j = 0; j < nyChunk; j++ )
                {
                    for ( int i = 0; i < nxChunk; i++ )
                    {
                        cout << "  " << P( chunkId, 1, i, j, k, 0 ) << '\t';
                        P( chunkId, 1, i, j, k, 0 ) = tmp[2 * ( nxChunk * nyChunk * k + nxChunk * j + i )];
                        P( chunkId, 1, i, j, k, 1 ) = tmp[2 * ( nxChunk * nyChunk * k + nxChunk * j + i ) + 1];
                    }
                    cout << endl;
                }
            }

            cout << "=============================" << endl;
            cout << "=============================" << endl;

            for ( int k = 0; k < nzChunk; k++ )
            {
    #if(OPENACC)
    #pragma acc loop vector collapse( 2 )
    #endif

                cout << "=============================" << endl;
                for ( int j = 0; j < nyChunk; j++ )
                {
                    for ( int i = 0; i < nxChunk; i++ )
                    {
                        cout << " " << tmp[2 * ( nxChunk * nyChunk * k + nxChunk * j + i )] << '\t';
                    }
                    cout << endl;
                }
            }
        }
    */
    switch ( next )
    {
        case 1:

            // rearranges x to y
            for ( int k = 0; k < nzChunk; k++ )
            {
                for ( int j = 0; j < nyChunk; j++ )
                {
                    for ( int i = 0; i < nxChunk; i++ )
                    {
                        P( chunkId, 1, i, j, k, 0 ) = tmp[2 * ( nxChunk * nyChunk * k + nxChunk * j + i )];
                        P( chunkId, 1, i, j, k, 1 ) = tmp[2 * ( nxChunk * nyChunk * k + nxChunk * j + i ) + 1];
                    }
                }
            }
            break;

        case 2:
            // rearrange y to z
            for ( int i = 0; i < nxChunk; i++ )
            {
                for ( int k = 0; k < nzChunk; k++ )
                {
                    for ( int j = 0; j < nyChunk; j++ )
                    {
                        P( chunkId, 2, i, j, k, 0 ) = tmp[2 * ( nxChunk * nyChunk * k + nyChunk * i + j )];
                        P( chunkId, 2, i, j, k, 1 ) = tmp[2 * ( nxChunk * nyChunk * k + nyChunk * i + j ) + 1];
                    }
                }
            }

            break;

        case 0:
            // rearrange z to x
            for ( int i = 0; i < nxChunk; i++ )
            {
                for ( int j = 0; j < nyChunk; j++ )
                {
                    for ( int k = 0; k < nzChunk; k++ )
                    {
                        P( chunkId, 0, i, j, k, 0 ) = tmp[2 * ( nzChunk * nyChunk * i + nzChunk * j + k )];
                        P( chunkId, 0, i, j, k, 1 ) = tmp[2 * ( nzChunk * nyChunk * i + nzChunk * j + k ) + 1];
                    }
                }
            }

            break;
    }
}
#endif

#if ( 0 )
#if ( OPENACC )
#pragma acc routine
#endif
void PencilDcmp::assignTemp( const int chunkId, const int original, const int next, double *tmp )
{

    switch ( next )
    {
        case 1:
            for ( int k = 0; k < nzChunk; k++ )
            {
#if ( OPENACC )
#pragma acc loop vector collapse( 2 )
#endif

                //        cout << "=============================" << endl;
                for ( int j = 0; j < nyChunk; j++ )
                {
                    for ( int i = 0; i < nxChunk; i++ )
                    {
                        tmp[2 * ( nxChunk * nyChunk * k + nxChunk * j + i )] = P( chunkId, original, i, j, k, 0 );
                        tmp[2 * ( nxChunk * nyChunk * k + nxChunk * j + i ) + 1] = P( chunkId, original, i, j, k, 1 );

                        //                cout << " " << tmp[2 * ( nxChunk * nyChunk * k + nxChunk * j + i )] << '\t';
                    }
                    //          cout << endl;
                }
            }
            break;
    }

    /*
        for ( int i = 0; i < 2 * nxChunk * nyChunk * nzChunk; i++ )
        {
            tmp[i] = P( chunkId * chunkSize + i );
            //     cout<<" myRank "<<myRank<< " " <<tmp[i] <<endl;
        }
    */
}
#endif

// new version does not use rearrange
#if ( 0 )
#if ( OPENACC )
#pragma acc routine gang
#endif
void PencilDcmp::rearrange( const int original, const int next )
{

#if ( OPENACC )
    double tmp[CHUNKSIZE];
#endif

    switch ( next )
    {
        case 1:

// rearranges x to y
#if ( OPENACC )
#pragma acc loop gang private( tmp[CHUNKSIZE] )
#endif
            for ( int id = 0; id < nChunk; id++ )
            {
                assignTemp( id, original, next, tmp );
#if ( OPENACC )
#pragma acc loop worker
#endif
                for ( int k = 0; k < nzChunk; k++ )
                {
#if ( OPENACC )
#pragma acc loop vector collapse( 2 )
#endif
                    for ( int j = 0; j < nyChunk; j++ )
                    {
                        for ( int i = 0; i < nxChunk; i++ )
                        {
                            P( id, 1, i, j, k, 0 ) = tmp[2 * ( nxChunk * nyChunk * k + nxChunk * j + i )];
                            P( id, 1, i, j, k, 1 ) = tmp[2 * ( nxChunk * nyChunk * k + nxChunk * j + i ) + 1];
                        }
                    }
                }
            }

            break;

        case 2:
// rearrange y to z
#if ( OPENACC )
#pragma acc loop gang private( tmp[CHUNKSIZE] )
#endif
            for ( int id = 0; id < nChunk; id++ )
            {
                assignTemp( id, original, next, tmp );
#if ( OPENACC )
#pragma acc loop worker
#endif
                for ( int i = 0; i < nxChunk; i++ )
                {
#if ( OPENACC )
#pragma acc loop vector collapse( 2 )
#endif
                    for ( int k = 0; k < nzChunk; k++ )
                    {
                        for ( int j = 0; j < nyChunk; j++ )
                        {
                            // P( id, 2, i, j, k, 0 ) = tmp[2 * ( nxChunk * nyChunk * k + nyChunk * i + j )];
                            P( id, 2, i, j, k, 0 ) = tmp[2 * ( nxChunk * nyChunk * k + nyChunk * i + j )];
                            P( id, 2, i, j, k, 1 ) = tmp[2 * ( nxChunk * nyChunk * k + nyChunk * i + j ) + 1];
                        }
                    }
                }
            }
            break;

        case 0:
// rearrange z to x
#if ( OPENACC )
#pragma acc loop gang private( tmp[CHUNKSIZE] )
#endif
            for ( int id = 0; id < nChunk; id++ )
            {
                assignTemp( id, original, next, tmp );
#if ( OPENACC )
#pragma acc loop worker
#endif
                for ( int i = 0; i < nxChunk; i++ )
                {
#if ( OPENACC )
#pragma acc loop vector collapse( 2 )
#endif
                    for ( int j = 0; j < nyChunk; j++ )
                    {
                        for ( int k = 0; k < nzChunk; k++ )
                        {
                            P( id, 0, i, j, k, 0 ) = tmp[2 * ( nzChunk * nyChunk * i + nzChunk * j + k )];
                            P( id, 0, i, j, k, 1 ) = tmp[2 * ( nzChunk * nyChunk * i + nzChunk * j + k ) + 1];
                        }
                    }
                }
            }

            break;
    }
}
#endif



    // xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

    //            Unecessary deleted methods, moved to junk

    // xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

    // void changeOwnership(); moved to junk.c

    //    void printX(); moved to junk.c

    //    void printY();
    //    sint newIndex(const sint chunkId,const sint j,const sint k);
    //    void decomposeNewIndex(const sint id, sint &chunkId, sint &j, sint &k);
    //    void rearrange( const int chunkId, const int original, const int next );
    //    void rearrange( const int original, const int next );
    //    void assignTemp( const int chunkId, const int original, const int next, double *tmp );
    //    void rearrangeY2Z(); /*!< the new algorithm does not use this as it performs Tridiagonal in-place */
    //    void assignTempY2Z( const int chunkId, const int k, double *tmp );
    //  void solver();

    //    void solverGPU();

    // void changeLocationXAux(sint i);
    // void initialize(); moved to junk.c

/*
#pragma acc routine gang
void Fourier::preprocessSignalDST00( const int direction )
{

#pragma acc loop gang
    for ( int k = 0; k < N; k++ )
    {
#pragma acc loop worker
        for ( int j = 0; j < N; j++ )
        {
            if ( direction == 0 )
            {
#pragma acc loop vector
                for ( int i = 0; i < N; i++ )
                {
                    Pn( N + i + 2, j, k, 1 ) = -Sn( N - 1 - i, j, k );
                    Pn( i + 1, j, k, 1 ) = Sn( i, j, k );
                }
            }
            else if ( direction == 1 )
            {
#pragma acc loop vector
                for ( int i = 0; i < N; i++ )
                {
                    Pn( j, N + i + 2, k, 1 ) = -Sn( j, N - 1 - i, k );
                    Pn( j, i + 1, k, 1 ) = Sn( j, i, k );
                }
            }
        }
    }
}
*/

#if(0)
void PencilDcmp::changeOwnershipPairwiseExchangeZX()
{

    MPI_Request request0, request1;
    MPI_Status status;

#if ( DEBUG_COMM )
    std::string filename = "rank";
    filename.append( to_string( myRank ) );
    ofstream myfile;
    myfile.open( filename );

#endif

#if ( DEBUG )
    int dest[2] = {myRank + 1, myRank - 1};
    for ( int i = 0; i < 2; i++ )
    {
        cout << " myRank  " << myRank << " will send to  " << getPeriodicRank( dest[i] ) << endl;
    }
#endif
    int i = myRank;
#if ( DEBUG )
    cout << " loop= " << ( p0 + 1 ) / 2 << endl;
#endif


    ///    P.moveHostToDevice( myRank );

  double *ptr1 = P.P;
  double *ptr2 = R.P;
#pragma acc host_data use_device( ptr1, ptr2 )
{
    for ( int j = 1; j < ( p0 ) / 2 + 1; j++ )
    {
       /// P.moveHostToDevice( getPeriodicIndex( i + j ) );

       ptr1= P.P+ P.chunkSize * getPeriodicIndex( i + j );
       MPI_Isend( ptr1 + P.chunkSize * getPeriodicIndex( i + j ), P.chunkSize, MPI_DOUBLE, getPeriodicRank( i + j ), myRank, Comm, &request0 );

       MPI_Irecv( ptr2, P.chunkSize, MPI_DOUBLE, getPeriodicRank( i - j ), getPeriodicRank( i - j ), Comm, &request1 );


        MPI_Wait( &request0, &status );

        /// once Isend completes move the chunk with (i-j) argument
        /// P.moveHostToDevice( getPeriodicIndex( i - j ) );
        /// switch send/recv order

        MPI_Isend( ptr1 + P.chunkSize * getPeriodicIndex( i - j ), P.chunkSize, MPI_DOUBLE, getPeriodicRank( i - j ), myRank, Comm,
                   &request0 );
        MPI_Irecv( ptr1+ P.chunkSize * getPeriodicIndex( i + j ), P.chunkSize, MPI_DOUBLE, getPeriodicRank( i + j ),
                   getPeriodicRank( i + j ), Comm, &request1 );

   ///     R.moveHostToDevice( );

        MPI_Wait( &request0, &status );
        MPI_Wait( &request1, &status );

#if ( DEBUG_COMM )

        myfile << "**********************************" << endl;
        myfile << " ZX myRank " << myRank << " destination " << getPeriodicRank( i - j ) << " tag and chunk id "
               << getPeriodicIndex( i - j ) << endl;
        myfile << " ZX myRank " << myRank << " source " << getPeriodicRank( i + j ) << " chunk ID " << getPeriodicIndex( i + j ) << endl;
        myfile << "**********************************" << endl;
#endif

        MPI_Wait( &request1, &status );
        for ( int k = 0; k < P.chunkSize; k++ )
        {
           P( k + P.chunkSize *getPeriodicIndex( i - j ) ) = R( k );
        }
    }
}
#if ( DEBUG )

    // cout<<"SUCCESS " <<endl;
    for ( int i = 0; i < nChunk; i++ )
    {
        cout << "ZX my rank " << myRank << " " << P( 0 + P.chunkSize * i ) << endl;
    }
#endif

#if ( DEBUG_COMM )
    myfile.close();
#endif
}
#endif


void PencilDcmp::nbrAllToAllYZ() /*!Two different communicators are required due to*/
{
// for All to allV
#if ( DEBUG )
    int *sndCnts = new int[p0];
    int *disp = new int[p0];

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

    // MPI_Neighbor_alltoallv(&P(0), sndCnts, int sdispls[], MPI_DOUBLE, recvbuf,
    // sndCnts, int rdispls[], MPI_DOUBLE ,nbrComm0, MPI_Request *request);

    MPI_Neighbor_alltoall( &P( 0 ), P.chunkSize, MPI_DOUBLE, &R( 0 ), P.chunkSize, MPI_DOUBLE, nbrComm0 );

    for ( int i = 0; i < nChunk * P.chunkSize; i++ )
    {
        P( i ) = R( i );
    }

    cout << " Rsize " << R.size() << " R chunk size " << R.chunkSize << endl;

#if ( DEBUG )
    for ( int i = 0; i < nChunk; i++ )
    {
        cout << "my rank " << myRank << " after transformation " << P( 0 + P.chunkSize * i ) << endl;
    }
#endif
}


#if ( 0 )
#if ( OPENACC )
#pragma acc routine gang
#endif
void PencilDcmp::rearrangeY2Z()
{
#if ( OPENACC )
    double tmp[CHUNKSIZE];
#endif
// rearranges x to y
#if ( OPENACC )
#pragma acc loop gang private( tmp[NZCHUNK *NXCHUNK] )
#endif
    for ( int id = 0; id < nChunk; id++ )
    {
#if ( OPENACC )
#pragma acc loop worker
#endif
        for ( int j = 0; j < nyChunk; j++ )
        {
            assignTempY2Z( id, j, tmp );
#if ( OPENACC )
#pragma acc loop vector collapse( 2 )
#endif
            for ( int k = 0; k < nzChunk; k++ )
            {
                for ( int i = 0; i < nxChunk; i++ )
                {
                    P( id, 1, i, j, k, 0 ) = tmp[2 * ( nxChunk * k + i )];
                    P( id, 1, i, j, k, 1 ) = tmp[2 * ( nxChunk * k + i ) + 1];
                }
            }
        }
    }
}

#if ( OPENACC )
#pragma acc routine worker
#endif
void PencilDcmp::assignTempY2Z( const int chunkId, const int j, double *tmp )
{

#if ( OPENACC )
#pragma acc loop worker
#endif
    for ( int k = 0; k < nzChunk; k++ )
    {

#if ( OPENACC )
#pragma acc loop vector
#endif
        for ( int i = 0; i < nxChunk; i++ )
        {
            tmp[2 * ( nxChunk * k + i )] = P( chunkId, 0, i, j, k, 0 );
            tmp[2 * ( nxChunk * k + i ) + 1] = P( chunkId, 0, i, j, k, 1 );
        }
    }
}
#endif

/*
inline void PencilDcmp::modifyEigForDirichlet( int i, int j, double *a )
{

    if ( j == 0 && faceTag[0] == 1 )
    {
        a[0] = -1.0;
    }
    if ( j == nxChunk - 1 && faceTag[1] == 1 )
    {
        a[0] = -1.0;
    }

    if ( i == 0 && faceTag[2] == 1 )
    {
        a[1] = -1.0;
    }
    if ( i == nyChunk - 1 && faceTag[3] == 1 )
    {
        a[1] = -1.0;
    }

   //  cout << " myRank "<<myRank<<" a[3]= "<< a[0] <<" "<<a[1] <<endl;
}
*/
/*
#elif( COMM_PATTERN == 0 )
void PencilDcmp::changeOwnershipPairwiseExchangeZX()
{

    MPI_Request request0, request1, request2, request3;
    MPI_Status status;

#if ( DEBUG_COMM )
    std::string filename = "rank";
    filename.append( to_string( myRank ) );
    ofstream myfile;
    myfile.open( filename );

#endif

#if ( DEBUG )
    int dest[2] = {myRank + 1, myRank - 1};
    for ( int i = 0; i < 2; i++ )
    {
        cout << " myRank  " << myRank << " will send to  " << getPeriodicRank( dest[i] ) << endl;
    }
#endif
    int i = myRank;
#if ( DEBUG )
    cout << " loop= " << ( p0 + 1 ) / 2 << endl;
#endif

    ///    P.moveHostToDevice( myRank );

    //  double *ptr1 = P.P;
    //  double *ptr2 = R.P;
    // #pragma acc host_data use_device( ptr1, ptr2 )
    {
        for ( int j = 1; j < ( p0 ) / 2 + 1; j++ )
        {
            /// P.moveHostToDevice( getPeriodicIndex( i + j ) );

            MPI_Irecv( &R( 0 ), P.chunkSize, MPI_DOUBLE, getPeriodicRank( i - j ), getPeriodicRank( i - j ), Comm, &request1 );

            MPI_Isend( &P( 0 ) + P.chunkSize * getPeriodicIndex( i + j ), P.chunkSize, MPI_DOUBLE, getPeriodicRank( i + j ), myRank, Comm,
                       &request0 );

            MPI_Wait( &request0, &status );

            /// once Isend completes move the chunk with (i-j) argument
            /// P.moveHostToDevice( getPeriodicIndex( i - j ) );
            /// switch send/recv order
            MPI_Irecv( &P( 0 ) + P.chunkSize * getPeriodicIndex( i + j ), P.chunkSize, MPI_DOUBLE, getPeriodicRank( i + j ),
                       getPeriodicRank( i + j ), Comm, &request3 );

            MPI_Isend( &P( 0 ) + P.chunkSize * getPeriodicIndex( i - j ), P.chunkSize, MPI_DOUBLE, getPeriodicRank( i - j ), myRank, Comm,
                       &request2 );
            ///     R.moveHostToDevice( );

            MPI_Wait( &request2, &status );
            MPI_Wait( &request3, &status );

#if ( DEBUG_COMM )

            myfile << "**********************************" << endl;
            myfile << " ZX myRank " << myRank << " destination " << getPeriodicRank( i - j ) << " tag and chunk id "
                   << getPeriodicIndex( i - j ) << endl;
            myfile << " ZX myRank " << myRank << " source " << getPeriodicRank( i + j ) << " chunk ID " << getPeriodicIndex( i + j )
                   << endl;
            myfile << "**********************************" << endl;
#endif

            MPI_Wait( &request1, &status );
            for ( int k = 0; k < P.chunkSize; k++ )
            {
                P( k + P.chunkSize *getPeriodicIndex( i - j ) ) = R( k );
            }
        }
    }
#if ( DEBUG )

    // cout<<"SUCCESS " <<endl;
    for ( int i = 0; i < nChunk; i++ )
    {
        cout << "ZX my rank " << myRank << " " << P( 0 + P.chunkSize * i ) << endl;
    }
#endif

#if ( DEBUG_COMM )
    myfile.close();
#endif
}
#endif
*/
/*
#if ( GPUAWARE2 == 0 && COMM_PATTERN == 0 )
void PencilDcmp::changeOwnershipPairwiseExchangeXY()
{
    // MPI_Request request[p0], request1[p0], request0[p0];
    MPI_Request request0, request1, request2, request3;
    MPI_Status status;
    int stride = p0;

// uses too much memory, use pairwise exhchange instead
//
// network congestion

#if ( DEBUG_COMM )
    std::string filename = "rank";
    filename.append( to_string( myRank ) );
    ofstream myfile;
    myfile.open( filename );

#endif

    int counter = 0;
    // MPI_Isend( &P(0), P.chunkSize, MPI_DOUBLE, getDestinationPeriodic(myRank),
    // myRank, Comm, &request1[0] );
    int i = myRank;
#if ( DEBUG )
    cout << " loop= " << ( p0 + 1 ) / 2 << endl;
#endif

    //  double *ptr1 = P.P;

    //  double *ptr2 = R.P;

    for ( int j = 1; j < ( p0 ) / 2 + 1; j++ )
    {

        MPI_Irecv( &R( 0 ), P.chunkSize, MPI_DOUBLE, getPeriodicRankStride( -j ), getPeriodicRankStride( -j ), Comm, &request1 );
        MPI_Isend( &P( 0 ) + P.chunkSize * getPeriodicIndexStride( j ), P.chunkSize, MPI_DOUBLE, getPeriodicRankStride( j ), myRank, Comm,
                   &request0 );
        MPI_Wait( &request0, &status );

#if ( DEBUG_COMM )
        myfile << "**********************************" << endl;
        myfile << " j " << j << " myRank " << myRank << " destination " << getPeriodicRankStride( j ) << " chunk id "
               << getPeriodicIndexStride( j ) << " tag " << myRank << endl;
        myfile << "**********************************" << endl;
        myfile << " j " << j << " myRank " << myRank << " source " << getPeriodicRankStride( -j ) << " chunk id "
               << getPeriodicIndexStride( -j ) << " tag " << getPeriodicRankStride( -j ) << endl;
#endif
        // write directly to the chunk that was previously send

        counter = 0;

        MPI_Irecv( &P( 0 ) + P.chunkSize * getPeriodicIndexStride( j ), P.chunkSize, MPI_DOUBLE, getPeriodicRankStride( j ),
                   getPeriodicRankStride( j ), Comm, &request3 );

        MPI_Isend( &P( 0 ) + P.chunkSize * getPeriodicIndexStride( -j ), P.chunkSize, MPI_DOUBLE, getPeriodicRankStride( -j ), myRank, Comm,
                   &request2 );

        MPI_Wait( &request2, &status );
        MPI_Wait( &request3, &status );

        // transfer from container to the desired location
        // cout<<nChunk<<endl;

        MPI_Wait( &request1, &status );

        for ( int k = 0; k < P.chunkSize; k++ )
        {
            P( k + P.chunkSize *getPeriodicIndexStride( -j ) ) = R( k );
        }

#if ( DEBUG_COMM )
        cout << " myRank " << myRank << " i " << i << " source " << getPeriodicRankStride( -j ) << " chunk ID "
             << getPeriodicIndexStride( -j ) << endl;
        cout << "**********************************" << endl;

        cout << " myRank " << myRank << " i " << i << " source " << getPeriodicRankStride( j ) << " chunk ID "
             << getPeriodicIndexStride( j ) << endl;
        cout << "**********************************" << endl;
#endif
    }

#if ( DEBUG_COMM )
    myfile.close();
#endif
}
*/
