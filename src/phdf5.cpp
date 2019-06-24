#include "phdf5.hpp"
#include "definitions.h"
//#include "stdio.h"
#define H5FILE_NAME "soln/Pxdmf3d%u.h5"
#define XDMF_NAME "soln/Pxdmf3d%u.xmf"
#define H5FILE "Pxdmf3d%u.h5"
#define DEBUG 0

void Phdf5::writeMultiBlock( ChunkedArray &F, uint appx )
{
    hid_t   file_id, dset_id;    /* file and dataset identifiers */
    hid_t   filespace, memspace; /* file and memory dataspace identifiers */
    hsize_t dimsf[4];            /* dataset dimensions */
    hsize_t chunk_dims[4];       /* chunk dimensions */
    hsize_t count[4];            /* hyperslab selection parameters */
    hsize_t block[4];
    hsize_t offset[4];
    hid_t   plist_id; /* property list identifier */
                      // uint    i, j, k, l;
    herr_t status;
    // int         *data=NULL;

    /*
     * MPI variables
     */
    // int           mpi_size, mpi_rank;
    PittPackReal *xtemp = NULL;
    PittPackReal *ytemp = NULL;
    PittPackReal *ztemp = NULL;
    PittPackReal *qtemp = NULL;

    uint L1 = F.nx;
    uint M1 = F.ny;
    uint N1 = F.nz;

    xtemp = new PittPackReal[L1 * M1 * N1];
    ytemp = new PittPackReal[L1 * M1 * N1];
    ztemp = new PittPackReal[L1 * M1 * N1];
    qtemp = new PittPackReal[L1 * M1 * N1];
    // cout<<L1<<"\t"<<M1<<"\t"<<N1<<endl;

    PittPackReal Xa;
    PittPackReal Xb;
    PittPackReal Ya;
    PittPackReal Yb;
    PittPackReal Za;
    PittPackReal Zb;

    char str[50];

    char str0[50];

    sprintf( str, H5FILE_NAME, appx );

    uint partialforestsize = F.nChunk; /*!<the forest size for each processor */

    CommPoint2Point<uint> com( &partialforestsize, 1 );
    uint                  offset1 = 0, totalvalue = 0;
    com.getOffset( partialforestsize, &offset1 );

    CommCollective<uint> comc( nullptr, 1 );
    comc.getTotalNumber( &offset1, &partialforestsize, &totalvalue );

    cout << " partial size " << partialforestsize << endl;

    hsize_t myoffset = offset1;
    //  total_size=totalvalue;

    totalnumber = totalvalue;

    cout << YELLOW "offsettttttttttttttt " << offset1 << RESET << endl;
    cout << YELLOW "L,M N " << L1 << " " << M1 << " " << N1 << RESET << endl;
    cout << YELLOW "totalvalue " << totalvalue << RESET << endl;

    /*
     * Set up file access property list with parallel I/O access
     */

    plist_id = H5Pcreate( H5P_FILE_ACCESS );
    H5Pset_fapl_mpio( plist_id, MPI_COMM_WORLD, MPI_INFO_NULL );

    /*
     * Create a new file collectively and release property list identifier.
     */

    file_id = H5Fcreate( str, H5F_ACC_TRUNC, H5P_DEFAULT, plist_id );
    H5Pclose( plist_id );

    /*
     * Create the dataspace for the dataset.
     */
    dimsf[0] = L1;
    dimsf[1] = M1;
    dimsf[2] = N1;
    dimsf[3] = totalnumber;

    chunk_dims[0] = L1;
    chunk_dims[1] = M1;
    chunk_dims[2] = N1;
    chunk_dims[3] = 1;

    /*
     * Each process defines dataset in memory and writes it to the hyperslab
     * in the file.
     */

    count[0] = 1;
    count[1] = 1;
    count[2] = 1;
    count[3] = 1;

    block[0] = chunk_dims[0];
    block[1] = chunk_dims[1];
    block[2] = chunk_dims[2];
    block[3] = chunk_dims[3];

    /*
      for loop such that each processor can write all the blocks
    */

    offset[0] = 0;
    offset[1] = 0;
    offset[2] = 0;

    // offset[3]=0/* will be set at loop*/
    //================================================================================
    //
    //                               write X
    //
    //================================================================================
    // the first argument is the dimension which is 4

    filespace = H5Screate_simple( 4, dimsf, NULL );
    memspace  = H5Screate_simple( 4, chunk_dims, NULL );

    sprintf( str0, "/X" );
    plist_id = H5Pcreate( H5P_DATASET_CREATE );
    H5Pset_chunk( plist_id, 4, chunk_dims );
    dset_id = H5Dcreate( file_id, str0, H5T_NATIVE_DOUBLE, filespace, H5P_DEFAULT, plist_id, H5P_DEFAULT );

    H5Pclose( plist_id );
    H5Sclose( filespace );

    filespace = H5Dget_space( dset_id );

    unsigned int level, coord_index;

    PittPackReal dx, dy, dz;
    PittPackReal denum;
    PittPackReal XYZ[6];
    //  morton<N> key;
    hsize_t co = 0;

    PittPackReal Xh, Yh, Zh;
    PittPackReal hx, hy, hz;
    int          index;
    /*
        cout << YELLOW << "Xa " << F.Xa << " Xb " << F.Xb << RESET << endl;
        cout << YELLOW << "Ya " << F.Ya << " Yb " << F.Yb << RESET << endl;
        cout << YELLOW << "Za " << F.Za << " Zb " << F.Zb << RESET << endl;
        cout << YELLOW << "Xh " << ( F.Xb - F.Xa ) / double( L1 - 1. ) << RESET << endl;
        cout << YELLOW << "F.nChunk " << F.nChunk << RESET << endl;
    */
    for ( int it = 0; it < F.nChunk; it++ )
    {
        Xa = F.Xa;
        Xb = F.Xb;
        hx = L1 - 1.0;
        Xh = ( Xb - Xa ) / ( hx );

        index = 0;

        for ( uint j = 0; j < L1; j++ )
        {
            for ( uint k = 0; k < M1; k++ )
            {
                for ( uint l = 0; l < N1; l++ )
                {
                    xtemp[index] = Xa + Xh * j;
                    index++;
                }
            }
        }
        // define the offset, only in the fourth dimension

        offset[3] = myoffset + co;

        printf( "offset=%llu\n", offset[3] );
        // status = H5Sselect_hyperslab(filespace, H5S_SELECT_SET, offset, NULL,
        // count, block);
        H5Sselect_hyperslab( filespace, H5S_SELECT_SET, offset, NULL, count, block );
        plist_id = H5Pcreate( H5P_DATASET_XFER );
        // H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);
        H5Pset_dxpl_mpio( plist_id, H5FD_MPIO_INDEPENDENT );
        // status = H5Dwrite(dset_id, H5T_NATIVE_DOUBLE,
        // memspace,filespace,plist_id, xtemp);

        if ( sizeof( PittPackReal ) == sizeof( double ) )
        {
            H5Dwrite( dset_id, H5T_NATIVE_DOUBLE, memspace, filespace, plist_id, xtemp );
        }
        else if ( sizeof( PittPackReal ) == sizeof( float ) )
        {
            H5Dwrite( dset_id, H5T_NATIVE_FLOAT, memspace, filespace, plist_id, xtemp );
        }

        H5Pclose( plist_id );
        co = co + 1;
    }

    H5Dclose( dset_id );

    // need to close dset_id since the new variable is going to be Y and Z ...
    // no need to open and close memspace and filespace all the time, simply reuse
    // them
    //================================================================================
    //
    //                               write Y
    //
    //================================================================================

    // filespace = H5Screate_simple(4, dimsf, NULL);
    //  memspace  = H5Screate_simple(4, chunk_dims, NULL);

    sprintf( str0, "/Y" );
    plist_id = H5Pcreate( H5P_DATASET_CREATE );
    H5Pset_chunk( plist_id, 4, chunk_dims );
    dset_id = H5Dcreate( file_id, str0, H5T_NATIVE_DOUBLE, filespace, H5P_DEFAULT, plist_id, H5P_DEFAULT );

    H5Pclose( plist_id );
    H5Sclose( filespace );

    filespace = H5Dget_space( dset_id );

    co = 0;

    for ( auto it = 0; it < F.nChunk; it++ )
    {
        //       proc.enclosingBox(key,XYZ);
        //       cout<<XYZ[0]<<" "<<XYZ[1]<<" "<<XYZ[2]<<" "<<XYZ[3]<<" "<<XYZ[4]<<"
        // "<<XYZ[5]<<endl;
        Ya = F.Ya;
        Yb = F.Yb;
        hy = M1 - 1.0;
        Yh = ( Yb - Ya ) / ( hy );

        index = 0;

        for ( uint j = 0; j < L1; j++ )
        {
            for ( uint k = 0; k < M1; k++ )
            {
                for ( uint l = 0; l < N1; l++ )
                {
                    ytemp[index] = Ya + Yh * k;
                    index++;
                }
            }
        }
        // define the offset, only in the fourth dimension

        offset[3] = myoffset + co;
        // printf("my_rank=%d offset=%d\n",my_rank,offset[3]);
        // status = H5Sselect_hyperslab(filespace, H5S_SELECT_SET, offset, NULL,
        // count, block);
        H5Sselect_hyperslab( filespace, H5S_SELECT_SET, offset, NULL, count, block );
        plist_id = H5Pcreate( H5P_DATASET_XFER );
        // H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);
        // status = H5Dwrite(dset_id, H5T_NATIVE_DOUBLE,
        // memspace,filespace,plist_id, xtemp);
        if ( sizeof( PittPackReal ) == sizeof( double ) )
        {
            H5Dwrite( dset_id, H5T_NATIVE_DOUBLE, memspace, filespace, plist_id, ytemp );
        }
        else if ( sizeof( PittPackReal ) == sizeof( float ) )
        {
            H5Dwrite( dset_id, H5T_NATIVE_FLOAT, memspace, filespace, plist_id, ytemp );
        }

        H5Pclose( plist_id );
        co = co + 1;
    }

    //================================================================================
    //
    //                               write Z
    //
    //================================================================================

    // filespace = H5Screate_simple(4, dimsf, NULL);
    //  memspace  = H5Screate_simple(4, chunk_dims, NULL);

    sprintf( str0, "/Z" );
    plist_id = H5Pcreate( H5P_DATASET_CREATE );
    H5Pset_chunk( plist_id, 4, chunk_dims );
    dset_id = H5Dcreate( file_id, str0, H5T_NATIVE_DOUBLE, filespace, H5P_DEFAULT, plist_id, H5P_DEFAULT );

    H5Pclose( plist_id );
    H5Sclose( filespace );

    filespace = H5Dget_space( dset_id );

    co = 0;

    PittPackReal Dz = ( F.Zb - F.Za ) / F.nChunk;

    for ( auto it = 0; it < F.nChunk; it++ )
    {
        // proc.enclosingBox(key,XYZ);
        Za = F.Za + it * Dz;
        Zb = F.Za + Dz * ( it + 1 );
        //      cout<<"Dz "<<Dz<<" Za "<<Za <<" Zb  "<<Zb<<endl;

        hz = N1 - 1.0;
        Zh = ( Zb - Za ) / ( hz );

        index = 0;

        for ( uint j = 0; j < L1; j++ )
        {
            for ( uint k = 0; k < M1; k++ )
            {
                for ( uint l = 0; l < N1; l++ )
                {
                    ztemp[index] = Za + Zh * l;
                    index++;
                }
            }
        }
        // define the offset, only in the fourth dimension

        offset[3] = myoffset + co;
        // printf("my_rank=%d offset=%d\n",my_rank,offset[3]);
        // status = H5Sselect_hyperslab(filespace, H5S_SELECT_SET, offset, NULL,
        // count, block);
        H5Sselect_hyperslab( filespace, H5S_SELECT_SET, offset, NULL, count, block );
        plist_id = H5Pcreate( H5P_DATASET_XFER );
        // H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);
        H5Pset_dxpl_mpio( plist_id, H5FD_MPIO_INDEPENDENT );
        // status = H5Dwrite(dset_id, H5T_NATIVE_DOUBLE,
        // memspace,filespace,plist_id, xtemp);
        if ( sizeof( PittPackReal ) == sizeof( double ) )
        {
            H5Dwrite( dset_id, H5T_NATIVE_DOUBLE, memspace, filespace, plist_id, ztemp );
        }
        else if ( sizeof( PittPackReal ) == sizeof( float ) )
        {
            H5Dwrite( dset_id, H5T_NATIVE_FLOAT, memspace, filespace, plist_id, ztemp );
        }
        H5Pclose( plist_id );
        co = co + 1;
    }
    //================================================================================
    //
    //                               write Q
    //
    //================================================================================
    // filespace = H5Screate_simple(4, dimsf, NULL);
    //  memspace  = H5Screate_simple(4, chunk_dims, NULL);

    sprintf( str0, "/Q" );
    plist_id = H5Pcreate( H5P_DATASET_CREATE );
    H5Pset_chunk( plist_id, 4, chunk_dims );
    dset_id = H5Dcreate( file_id, str0, H5T_NATIVE_DOUBLE, filespace, H5P_DEFAULT, plist_id, H5P_DEFAULT );

    H5Pclose( plist_id );
    H5Sclose( filespace );

    filespace = H5Dget_space( dset_id );

    co = 0;

    int my_rank, comsize;
    MPI_Comm_rank( MPI_COMM_WORLD, &my_rank );
    MPI_Comm_size( MPI_COMM_WORLD, &comsize );

    for ( int it = 0; it < F.nChunk; it++ )
    {
        index = 0;

        for ( uint j = 0; j < L1; j++ )
        {
            for ( uint k = 0; k < M1; k++ )
            {
                for ( uint l = 0; l < N1; l++ )
                {
                    // ztemp[index]=com.myRank();
                    //                  ztemp[index]=F.p(2*index);
                    //                    qtemp[index] = my_rank;
                    qtemp[index] = F( it, 0, j, k, l, 0 );
                    index++;
                }
            }
        }
        // define the offset, only in the fourth dimension

        offset[3] = myoffset + co;
        // printf("my_rank=%d offset=%d\n",my_rank,offset[3]);
        // status = H5Sselect_hyperslab(filespace, H5S_SELECT_SET, offset, NULL,
        // count, block);
        H5Sselect_hyperslab( filespace, H5S_SELECT_SET, offset, NULL, count, block );
        plist_id = H5Pcreate( H5P_DATASET_XFER );
        // H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);
        H5Pset_dxpl_mpio( plist_id, H5FD_MPIO_INDEPENDENT );
        // status = H5Dwrite(dset_id, H5T_NATIVE_DOUBLE,
        // memspace,filespace,plist_id, xtemp);

        if ( sizeof( PittPackReal ) == sizeof( double ) )
        {
            H5Dwrite( dset_id, H5T_NATIVE_DOUBLE, memspace, filespace, plist_id, qtemp );
            //   cout<<BLUE "using double " RESET<<endl;
        }
        else if ( sizeof( PittPackReal ) == sizeof( float ) )
        {
            H5Dwrite( dset_id, H5T_NATIVE_FLOAT, memspace, filespace, plist_id, qtemp );
        }

        //      H5Dwrite(dset_id, H5T_NATIVE_DOUBLE, memspace,filespace,plist_id,
        // ztemp);
        H5Pclose( plist_id );
        co = co + 1;
    }

    xdmfMultiBlock( F, com.mySize(), com.myRank(), offset1, appx );

    delete[] xtemp;
    delete[] ytemp;
    delete[] ztemp;
    delete[] qtemp;
}

static void integer_string( char *strin, int i )
{
    if ( i < 10 )
    {
        sprintf( strin, "00000%d", i );
    }
    else if ( i < 100 )
    {
        sprintf( strin, "0000%d", i );
    }
    else if ( i < 1000 )
    {
        sprintf( strin, "000%d", i );
    }
    else if ( i < 10000 )
    {
        sprintf( strin, "00%d", i );
    }
    else if ( i < 100000 )
    {
        sprintf( strin, "0%d", i );
    }
    else if ( i < 1000000 )
    {
        sprintf( strin, "%d", i );
    }
    else
    {
        printf( "not in the range, offset is too big\n" );
        exit( 0 );
    }
}

//#define offset0 156
//#define offset1 2005

void Phdf5::xdmfMultiBlock( ChunkedArray &F, integer comsize, integer my_rank, uint offset, uint appx )
{
    MPI_File     fp;
    int          buf[1000], np = comsize;
    MPI_Request  request;
    unsigned int i;
    int          j;

    MPI_Status  status;
    const char *names[] = {"X", "Y", "Z"};

    uint L1 = F.nx;
    uint M1 = F.ny;
    uint N1 = F.nz;

    char strname[1000], fname[1000];
    uint ncube_total = totalnumber;
    sprintf( strname, XDMF_NAME, appx );
    char str[1000];

    const int offset0 = 156;
    const int offset1 = 2005;

    //   cout << "size str" << strlen( strname ) << endl;

    MPI_File_open( MPI_COMM_WORLD, strname, MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &fp );

    // printf("%s\n",str);
    // dont need to specify the address here as they are both reside in one
    // folder, that is
    // why I use H5FILE in addition to H5FILE_NAME
    sprintf( fname, H5FILE, appx );

    char strL[100];
    char strM[100];
    char strN[100];
    char strNcube[1000];
    char stroff[1000];
    int  index;

    if ( L1 < 10 )
    {
        sprintf( strL, "00%d", L1 );
        sprintf( strM, "00%d", M1 );
        sprintf( strN, "00%d", N1 );
    }
    else if ( L1 < 100 )
    {
        sprintf( strL, "0%d", L1 );
        sprintf( strM, "0%d", M1 );
        sprintf( strN, "0%d", N1 );
    }
    else if ( L1 < 1000 )
    {
        sprintf( strL, "%d", L1 );
        sprintf( strM, "%d", M1 );
        sprintf( strN, "%d", N1 );
    }

    else
    {
        printf( "discretization too big go change xdmf function\n" );
        exit( 0 );
    }

    if ( ncube_total < 10 )
    {
        sprintf( strNcube, "00000%d", ncube_total );
    }
    else if ( ncube_total < 100 )
    {
        sprintf( strNcube, "0000%d", ncube_total );
    }
    else if ( ncube_total < 1000 )
    {
        sprintf( strNcube, "000%d", ncube_total );
    }
    else if ( ncube_total < 10000 )
    {
        sprintf( strNcube, "00%d", ncube_total );
    }
    else if ( ncube_total < 100000 )
    {
        sprintf( strNcube, "0%d", ncube_total );
    }
    else if ( ncube_total < 1000000 )
    {
        sprintf( strNcube, "%d", ncube_total );
    }
    else
    {
        printf( "number of cubes larger than 1000000\n" );
        exit( 0 );
    }

    cout << ncube_total << endl;
    uint co = 0;
    // a counts the offset for header whic is only written by process rank 0 and
    // and b the hyperslab part for each cube
    int a = 0, b = 0;
    if ( my_rank == 0 )
    {
        sprintf( str, "<?xml version=\"1.0\" ?>\n" );
        MPI_File_write( fp, str, strlen( str ), MPI_CHAR, &status );
        a = a + strlen( str );
        sprintf( str, "<!DOCTYPE Xdmf SYSTEM \"Xdmf.dtd\" []> \n" );
        MPI_File_write( fp, str, strlen( str ), MPI_CHAR, &status );
        a = a + strlen( str );
        sprintf( str, "<Xdmf Version=\"2.0\">\n" );
        MPI_File_write( fp, str, strlen( str ), MPI_CHAR, &status );
        a = a + strlen( str );
        sprintf( str, "<Domain>\n" );
        MPI_File_write( fp, str, strlen( str ), MPI_CHAR, &status );
        a = a + strlen( str );
        sprintf( str, "<Grid Name=\"AMR\" GridType=\"Collection\" "
                      "CollectionType=\"Spatial\">\n" );
        MPI_File_write( fp, str, strlen( str ), MPI_CHAR, &status );
        a = a + strlen( str );
        printf( "size a=%d\n", a );
        // for(unsigned int i=0;i<cube.size();i++)
        co = 0;

        // for(auto it=F.trees.begin();it!=F.trees.end();it++)
        {
            //        for(auto it2=(*it).begin();it2!=(*it).end();it2++)
            for ( auto it2 = 0; it2 < F.nChunk; it2++ )
            {
                // b calculates the offset for each block, need to be set to zero
                // becasue it is inside the loop
                // but we dont want to accumulate, calculate only for one block
                b = 0;
                sprintf( str, "   <Grid Name=\"mesh0\" GridType=\"Uniform\">\n" );
                MPI_File_write( fp, str, strlen( str ), MPI_CHAR, &status );
                b = b + strlen( str );
                sprintf( str,
                         "        <Topology TopologyType=\"3DSMesh\" "
                         "NumberOfElements=\"%s %s %s\"/>\n",
                         strL, strM, strN );
                // printf("hyperslab =%d\n",strlen(str));
                MPI_File_write( fp, str, strlen( str ), MPI_CHAR, &status );
                b = b + strlen( str );
                sprintf( str, "          <Geometry GeometryType=\"X_Y_Z\">  \n" );
                MPI_File_write( fp, str, strlen( str ), MPI_CHAR, &status );
                b = b + strlen( str );

                for ( j = 0; j < 3; j++ )
                {
                    sprintf( str,
                             "      <DataItem ItemType=\"HyperSlab\" Dimensions=\"%s "
                             "%s %s %d\" NumberType=\"Float\" Precision=\"4\"  "
                             "Type=\"HyperSlab\"> \n",
                             strL, strM, strN, 1 );
                    MPI_File_write( fp, str, strlen( str ), MPI_CHAR, &status );
                    b = b + strlen( str );
                    sprintf( str, "         <DataItem Dimensions=\"3 4\" Format=\"XML\" > \n" );
                    MPI_File_write( fp, str, strlen( str ), MPI_CHAR, &status );
                    b     = b + strlen( str );
                    index = offset + co;
                    integer_string( stroff, index );
                    sprintf( str, "         %d %d %d %s  \n", 0, 0, 0, stroff );
                    MPI_File_write( fp, str, strlen( str ), MPI_CHAR, &status );
                    b = b + strlen( str );
                    sprintf( str, "         %d %d %d %d  \n", 1, 1, 1, 1 );
                    MPI_File_write( fp, str, strlen( str ), MPI_CHAR, &status );
                    b = b + strlen( str );
                    sprintf( str, "         %s %s %s %d  \n", strL, strM, strN, 1 );
                    MPI_File_write( fp, str, strlen( str ), MPI_CHAR, &status );
                    b = b + strlen( str );
                    sprintf( str, "         </DataItem>\n" );
                    MPI_File_write( fp, str, strlen( str ), MPI_CHAR, &status );
                    b = b + strlen( str );
                    sprintf( str,
                             "          <DataItem Name=\"%s\" Dimensions=\"%s %s %s "
                             "%s\" NumberType=\"Float\" Precision=\"4\" "
                             "Format=\"HDF\">\n",
                             names[j], strL, strM, strN, strNcube );
                    MPI_File_write( fp, str, strlen( str ), MPI_CHAR, &status );
                    b = b + strlen( str );
                    sprintf( str, "          %s:/%s\n", fname, names[j] );
                    MPI_File_write( fp, str, strlen( str ), MPI_CHAR, &status );
                    b = b + strlen( str );
                    sprintf( str, "         </DataItem>\n" );
                    MPI_File_write( fp, str, strlen( str ), MPI_CHAR, &status );
                    b = b + strlen( str );
                    sprintf( str, "         </DataItem>\n" );
                    MPI_File_write( fp, str, strlen( str ), MPI_CHAR, &status );
                    b = b + strlen( str );
                }
                sprintf( str, "      </Geometry>   \n" );
                MPI_File_write( fp, str, strlen( str ), MPI_CHAR, &status );
                b = b + strlen( str );
                // added because of including solution vector Q
                sprintf( str, "         <Attribute Name=\"Q\" AttributeType=\"Scalar\" "
                              "Center=\"Node\"> \n" );
                MPI_File_write( fp, str, strlen( str ), MPI_CHAR, &status );
                b = b + strlen( str );
                sprintf( str,
                         "      <DataItem ItemType=\"HyperSlab\" Dimensions=\"%s "
                         "%s %s %d\" NumberType=\"Float\" Precision=\"4\"  "
                         "Type=\"HyperSlab\"> \n",
                         strL, strM, strN, 1 );
                MPI_File_write( fp, str, strlen( str ), MPI_CHAR, &status );
                b = b + strlen( str );
                sprintf( str, "         <DataItem Dimensions=\"3 4\" Format=\"XML\" > \n" );
                MPI_File_write( fp, str, strlen( str ), MPI_CHAR, &status );
                b     = b + strlen( str );
                index = offset + co;
                integer_string( stroff, index );
                sprintf( str, "         %d %d %d %s  \n", 0, 0, 0, stroff );
                MPI_File_write( fp, str, strlen( str ), MPI_CHAR, &status );
                b = b + strlen( str );
                sprintf( str, "         %d %d %d %d  \n", 1, 1, 1, 1 );
                MPI_File_write( fp, str, strlen( str ), MPI_CHAR, &status );
                b = b + strlen( str );
                sprintf( str, "         %s %s %s %d  \n", strL, strM, strN, 1 );
                MPI_File_write( fp, str, strlen( str ), MPI_CHAR, &status );
                b = b + strlen( str );
                sprintf( str, "         </DataItem>\n" );
                MPI_File_write( fp, str, strlen( str ), MPI_CHAR, &status );
                b = b + strlen( str );
                sprintf( str,
                         "          <DataItem Name=\"Q\" Dimensions=\"%s %s %s "
                         "%s\" NumberType=\"Float\" Precision=\"4\" "
                         "Format=\"HDF\">\n",
                         strL, strM, strN, strNcube );
                MPI_File_write( fp, str, strlen( str ), MPI_CHAR, &status );
                b = b + strlen( str );
                sprintf( str, "          %s:/Q\n", fname );
                MPI_File_write( fp, str, strlen( str ), MPI_CHAR, &status );
                b = b + strlen( str );

                sprintf( str, "         </DataItem>\n" );
                MPI_File_write( fp, str, strlen( str ), MPI_CHAR, &status );
                b = b + strlen( str );
                sprintf( str, "         </DataItem>\n" );
                MPI_File_write( fp, str, strlen( str ), MPI_CHAR, &status );
                b = b + strlen( str );
                sprintf( str, " </Attribute>\n" );
                MPI_File_write( fp, str, strlen( str ), MPI_CHAR, &status );
                b = b + strlen( str );
                sprintf( str, "  </Grid>\n" );
                MPI_File_write( fp, str, strlen( str ), MPI_CHAR, &status );
                b = b + strlen( str );
#if ( 1 )
                co++;
                // procesor zero should broadcast this value to others
                // since this is inside an if loop we can not use collective
                // rather we simply send via send and recieves
                // put error here
                if ( a != offset0 || b != offset1 )
                {
                    printf( "??????????????????????????????????????????\n" );
                    printf( "Go fix your offset for Xdmf meta data file a=%d offset0=%d "
                            "b=%d offset1=%d\n",
                            a, offset0, b, offset1 );
                    printf( "??????????????????????????????????????????\n" );
                    exit( 0 );
                }

//              printf("Go fix your offset for Xdmf meta data file a=%d
// offset0=%d b=%d offset1=%d\n",a,offset0,b,offset1);
#endif
            }
        }
    }
#if ( 1 )
    else
    {
        // recieve the offset value form proc 0 and set your file offse taccordingly
        // other ranks do something else
        // 5 is added because of the header of the file written by  processor 0

        int mpi_offset = 0;
        co             = 0;

        for ( int it = 0; it < F.nChunk; it++ )
        {
            mpi_offset = ( ( offset + co ) * offset1 + offset0 );
            mpi_offset = mpi_offset * sizeof( char );
            // printf("my_rank = %d offset =%d mpi_offset=%d
            // i=%d\n",my_rank,offset,mpi_offset,i);

            MPI_File_seek( fp, mpi_offset, MPI_SEEK_SET );
            // MPI_File_set_view(fp,mpi_offset,MPI_CHAR,MPI_CHAR, "native",
            // MPI_INFO_NULL );
            sprintf( str, "   <Grid Name=\"mesh0\" GridType=\"Uniform\">\n" );
            MPI_File_write( fp, str, strlen( str ), MPI_CHAR, &status );
            sprintf( str,
                     "        <Topology TopologyType=\"3DSMesh\" "
                     "NumberOfElements=\"%s %s %s\"/>\n",
                     strL, strM, strN );
            MPI_File_write( fp, str, strlen( str ), MPI_CHAR, &status );
            sprintf( str, "          <Geometry GeometryType=\"X_Y_Z\">  \n" );
            MPI_File_write( fp, str, strlen( str ), MPI_CHAR, &status );

#if ( 1 )
            for ( j = 0; j < 3; j++ )
            {
                sprintf( str,
                         "      <DataItem ItemType=\"HyperSlab\" Dimensions=\"%s "
                         "%s %s %d\" NumberType=\"Float\" Precision=\"4\"  "
                         "Type=\"HyperSlab\"> \n",
                         strL, strM, strN, 1 );
                MPI_File_write( fp, str, strlen( str ), MPI_CHAR, &status );
                sprintf( str, "         <DataItem Dimensions=\"3 4\" Format=\"XML\" > \n" );
                MPI_File_write( fp, str, strlen( str ), MPI_CHAR, &status );
                index = offset + co;
                integer_string( stroff, index );
                sprintf( str, "         %d %d %d %s  \n", 0, 0, 0, stroff );
                MPI_File_write( fp, str, strlen( str ), MPI_CHAR, &status );
                sprintf( str, "         %d %d %d %d  \n", 1, 1, 1, 1 );
                MPI_File_write( fp, str, strlen( str ), MPI_CHAR, &status );
                sprintf( str, "         %s %s %s %d  \n", strL, strM, strN, 1 );
                MPI_File_write( fp, str, strlen( str ), MPI_CHAR, &status );
                sprintf( str, "         </DataItem>\n" );
                MPI_File_write( fp, str, strlen( str ), MPI_CHAR, &status );
                sprintf( str,
                         "          <DataItem Name=\"%s\" Dimensions=\"%s %s %s "
                         "%s\" NumberType=\"Float\" Precision=\"4\" "
                         "Format=\"HDF\">\n",
                         names[j], strL, strM, strN, strNcube );
                MPI_File_write( fp, str, strlen( str ), MPI_CHAR, &status );
                sprintf( str, "          %s:/%s\n", fname, names[j] );
                MPI_File_write( fp, str, strlen( str ), MPI_CHAR, &status );
                sprintf( str, "         </DataItem>\n" );
                MPI_File_write( fp, str, strlen( str ), MPI_CHAR, &status );
                sprintf( str, "         </DataItem>\n" );
                MPI_File_write( fp, str, strlen( str ), MPI_CHAR, &status );
            }
#endif
            sprintf( str, "      </Geometry>   \n" );
            MPI_File_write( fp, str, strlen( str ), MPI_CHAR, &status );

// insert here the value of qs
#if ( 1 )
            sprintf( str, "         <Attribute Name=\"Q\" AttributeType=\"Scalar\" "
                          "Center=\"Node\"> \n" );
            MPI_File_write( fp, str, strlen( str ), MPI_CHAR, &status );
            sprintf( str,
                     "      <DataItem ItemType=\"HyperSlab\" Dimensions=\"%s %s "
                     "%s %d\" NumberType=\"Float\" Precision=\"4\"  "
                     "Type=\"HyperSlab\"> \n",
                     strL, strM, strN, 1 );
            MPI_File_write( fp, str, strlen( str ), MPI_CHAR, &status );
            sprintf( str, "         <DataItem Dimensions=\"3 4\" Format=\"XML\" > \n" );
            MPI_File_write( fp, str, strlen( str ), MPI_CHAR, &status );
            index = offset + co;
            integer_string( stroff, index );
            sprintf( str, "         %d %d %d %s  \n", 0, 0, 0, stroff );
            MPI_File_write( fp, str, strlen( str ), MPI_CHAR, &status );
            sprintf( str, "         %d %d %d %d  \n", 1, 1, 1, 1 );
            MPI_File_write( fp, str, strlen( str ), MPI_CHAR, &status );
            sprintf( str, "         %s %s %s %d  \n", strL, strM, strN, 1 );
            MPI_File_write( fp, str, strlen( str ), MPI_CHAR, &status );
            sprintf( str, "         </DataItem>\n" );
            MPI_File_write( fp, str, strlen( str ), MPI_CHAR, &status );
            sprintf( str,
                     "          <DataItem Name=\"Q\" Dimensions=\"%s %s %s %s\" "
                     "NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n",
                     strL, strM, strN, strNcube );
            MPI_File_write( fp, str, strlen( str ), MPI_CHAR, &status );
            sprintf( str, "          %s:/Q\n", fname );
            MPI_File_write( fp, str, strlen( str ), MPI_CHAR, &status );
            sprintf( str, "         </DataItem>\n" );
            MPI_File_write( fp, str, strlen( str ), MPI_CHAR, &status );
            sprintf( str, "         </DataItem>\n" );
            MPI_File_write( fp, str, strlen( str ), MPI_CHAR, &status );
            sprintf( str, " </Attribute>\n" );
            MPI_File_write( fp, str, strlen( str ), MPI_CHAR, &status );

            sprintf( str, "  </Grid>\n" );
            MPI_File_write( fp, str, strlen( str ), MPI_CHAR, &status );
            co++;
#endif
        }
        // add the last part of the xdmf text file to close the arguments
    }
#endif
#if ( 1 )
    if ( my_rank == ( np - 1 ) )
    {
        sprintf( str, "  </Grid>\n" );
        MPI_File_write( fp, str, strlen( str ), MPI_CHAR, &status );
        sprintf( str, "      </Domain>   \n" );
        MPI_File_write( fp, str, strlen( str ), MPI_CHAR, &status );
        sprintf( str, "  </Xdmf>\n" );
        MPI_File_write( fp, str, strlen( str ), MPI_CHAR, &status );
    }
#endif
    MPI_File_close( &fp );
}
#if ( DEBUG )
void Phdf5::writeMultiBlockCellCenter( ChunkedArray &F, uint appx )
{
    hid_t   file_id, dset_id;    /* file and dataset identifiers */
    hid_t   filespace, memspace; /* file and memory dataspace identifiers */
    hsize_t dimsf[4];            /* dataset dimensions */
    hsize_t chunk_dims[4];       /* chunk dimensions */
    hsize_t count[4];            /* hyperslab selection parameters */
    hsize_t block[4];
    hsize_t offset[4];
    hid_t   plist_id; /* property list identifier */
    uint    i, j, k, l;
    herr_t  status;
    // int         *data=NULL;

    /*
     * MPI variables
     */
    int           mpi_size, mpi_rank;
    PittPackReal *xtemp = NULL;
    PittPackReal *ytemp = NULL;
    PittPackReal *ztemp = NULL;
    PittPackReal *qtemp = NULL;

    uint L1 = F.nx + 1;
    uint M1 = F.ny + 1;
    uint N1 = F.nz + 1;

    xtemp = new PittPackReal[L1 * M1 * N1];
    ytemp = new PittPackReal[L1 * M1 * N1];
    ztemp = new PittPackReal[L1 * M1 * N1];
    qtemp = new PittPackReal[( L1 - 1 ) * ( M1 - 1 ) * ( N1 - 1 )];
    // cout<<L1<<"\t"<<M1<<"\t"<<N1<<endl;

    PittPackReal Xa;
    PittPackReal Xb;
    PittPackReal Ya;
    PittPackReal Yb;
    PittPackReal Za;
    PittPackReal Zb;

    char str[50];

    char str0[50];

    sprintf( str, H5FILE_NAME, appx );

    uint partialforestsize = F.nChunk; /*!<the forest size for each processor */

    CommPoint2Point<uint> com( &partialforestsize, 1 );
    uint                  offset1 = 0, totalvalue = 0;
    com.getOffset( partialforestsize, &offset1 );

    CommCollective<uint> comc( nullptr, 1 );
    comc.getTotalNumber( &offset1, &partialforestsize, &totalvalue );

    cout << " partial size " << partialforestsize << endl;

    hsize_t myoffset = offset1;
    //  total_size=totalvalue;

    totalnumber = totalvalue;
    /*
        cout << YELLOW "offsettttttttttttttt " << offset1 << RESET << endl;
        cout << YELLOW "L,M N " << L1 << " " << M1 << " " << N1 << RESET << endl;
        cout << YELLOW "totalvalue " << totalvalue << RESET << endl;
    */
    /*
     * Set up file access property list with parallel I/O access
     */

    plist_id = H5Pcreate( H5P_FILE_ACCESS );
    H5Pset_fapl_mpio( plist_id, MPI_COMM_WORLD, MPI_INFO_NULL );

    /*
     * Create a new file collectively and release property list identifier.
     */

    file_id = H5Fcreate( str, H5F_ACC_TRUNC, H5P_DEFAULT, plist_id );
    H5Pclose( plist_id );

    /*
     * Create the dataspace for the dataset.
     */
    dimsf[0] = L1;
    dimsf[1] = M1;
    dimsf[2] = N1;
    dimsf[3] = totalnumber;

    chunk_dims[0] = L1;
    chunk_dims[1] = M1;
    chunk_dims[2] = N1;
    chunk_dims[3] = 1;

    /*
     * Each process defines dataset in memory and writes it to the hyperslab
     * in the file.
     */

    count[0] = 1;
    count[1] = 1;
    count[2] = 1;
    count[3] = 1;

    block[0] = chunk_dims[0];
    block[1] = chunk_dims[1];
    block[2] = chunk_dims[2];
    block[3] = chunk_dims[3];

    /*
      for loop such that each processor can write all the blocks
    */
    offset[0] = 0;
    offset[1] = 0;
    offset[2] = 0;
    // offset[3] will be set at loop
    //================================================================================
    //
    //                               write X
    //
    //================================================================================
    // the first argument is the dimension which is 4

    filespace = H5Screate_simple( 4, dimsf, NULL );
    memspace  = H5Screate_simple( 4, chunk_dims, NULL );

    sprintf( str0, "/X" );
    plist_id = H5Pcreate( H5P_DATASET_CREATE );
    H5Pset_chunk( plist_id, 4, chunk_dims );
    dset_id = H5Dcreate( file_id, str0, H5T_NATIVE_DOUBLE, filespace, H5P_DEFAULT, plist_id, H5P_DEFAULT );

    H5Pclose( plist_id );
    H5Sclose( filespace );

    filespace = H5Dget_space( dset_id );

    unsigned int level, coord_index;

    PittPackReal dx, dy, dz;
    PittPackReal denum;
    PittPackReal XYZ[6];
    //  morton<N> key;
    hsize_t co = 0;

    PittPackReal Xh, Yh, Zh;
    PittPackReal hx, hy, hz;
    int          index;
    /*
        cout << YELLOW << "Xa " << F.Xa << " Xb " << F.Xb << RESET << endl;
        cout << YELLOW << "Ya " << F.Ya << " Yb " << F.Yb << RESET << endl;
        cout << YELLOW << "Za " << F.Za << " Zb " << F.Zb << RESET << endl;
        cout << YELLOW << "Xh " << ( F.Xb - F.Xa ) / double( L1 - 1. ) << RESET << endl;
        cout << YELLOW << "F.nChunk " << F.nChunk << RESET << endl;
    */
    for ( int it = 0; it < F.nChunk; it++ )
    {
        Xa = F.Xa;
        Xb = F.Xb;
        hx = L1 - 1.0;
        Xh = ( Xb - Xa ) / ( hx );

        index = 0;

        for ( uint j = 0; j < L1; j++ )
        {
            for ( uint k = 0; k < M1; k++ )
            {
                for ( uint l = 0; l < N1; l++ )
                {
                    xtemp[index] = Xa + Xh * j;
                    index++;
                }
            }
        }
        // define the offset, only in the fourth dimension

        offset[3] = myoffset + co;

        printf( "offset=%d\n", offset[3] );
        // status = H5Sselect_hyperslab(filespace, H5S_SELECT_SET, offset, NULL,
        // count, block);
        H5Sselect_hyperslab( filespace, H5S_SELECT_SET, offset, NULL, count, block );
        plist_id = H5Pcreate( H5P_DATASET_XFER );
        // H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);
        H5Pset_dxpl_mpio( plist_id, H5FD_MPIO_INDEPENDENT );
        // status = H5Dwrite(dset_id, H5T_NATIVE_DOUBLE,
        // memspace,filespace,plist_id, xtemp);

        if ( sizeof( PittPackReal ) == sizeof( double ) )
        {
            H5Dwrite( dset_id, H5T_NATIVE_DOUBLE, memspace, filespace, plist_id, xtemp );
        }
        else if ( sizeof( PittPackReal ) == sizeof( float ) )
        {
            H5Dwrite( dset_id, H5T_NATIVE_FLOAT, memspace, filespace, plist_id, xtemp );
        }

        H5Pclose( plist_id );
        co = co + 1;
    }

    H5Dclose( dset_id );

    // need to close dset_id since the new variable is going to be Y and Z ...
    // no need to open and close memspace and filespace all the time, simply reuse
    // them
    //================================================================================
    //
    //                               write Y
    //
    //================================================================================

    // filespace = H5Screate_simple(4, dimsf, NULL);
    //  memspace  = H5Screate_simple(4, chunk_dims, NULL);

    sprintf( str0, "/Y" );
    plist_id = H5Pcreate( H5P_DATASET_CREATE );
    H5Pset_chunk( plist_id, 4, chunk_dims );
    dset_id = H5Dcreate( file_id, str0, H5T_NATIVE_DOUBLE, filespace, H5P_DEFAULT, plist_id, H5P_DEFAULT );

    H5Pclose( plist_id );
    H5Sclose( filespace );

    filespace = H5Dget_space( dset_id );

    co = 0;

    for ( auto it = 0; it < F.nChunk; it++ )
    {
        //       proc.enclosingBox(key,XYZ);
        //       cout<<XYZ[0]<<" "<<XYZ[1]<<" "<<XYZ[2]<<" "<<XYZ[3]<<" "<<XYZ[4]<<"
        // "<<XYZ[5]<<endl;
        Ya = F.Ya;
        Yb = F.Yb;
        hy = M1 - 1.0;
        Yh = ( Yb - Ya ) / ( hy );

        index = 0;

        for ( uint j = 0; j < L1; j++ )
        {
            for ( uint k = 0; k < M1; k++ )
            {
                for ( uint l = 0; l < N1; l++ )
                {
                    ytemp[index] = Ya + Yh * k;
                    index++;
                }
            }
        }
        // define the offset, only in the fourth dimension

        offset[3] = myoffset + co;
        // printf("my_rank=%d offset=%d\n",my_rank,offset[3]);
        // status = H5Sselect_hyperslab(filespace, H5S_SELECT_SET, offset, NULL,
        // count, block);
        H5Sselect_hyperslab( filespace, H5S_SELECT_SET, offset, NULL, count, block );
        plist_id = H5Pcreate( H5P_DATASET_XFER );
        // H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);
        // status = H5Dwrite(dset_id, H5T_NATIVE_DOUBLE,
        // memspace,filespace,plist_id, xtemp);
        if ( sizeof( PittPackReal ) == sizeof( double ) )
        {
            H5Dwrite( dset_id, H5T_NATIVE_DOUBLE, memspace, filespace, plist_id, ytemp );
        }
        else if ( sizeof( PittPackReal ) == sizeof( float ) )
        {
            H5Dwrite( dset_id, H5T_NATIVE_FLOAT, memspace, filespace, plist_id, ytemp );
        }

        H5Pclose( plist_id );
        co = co + 1;
    }

    //================================================================================
    //
    //                               write Z
    //
    //================================================================================

    // filespace = H5Screate_simple(4, dimsf, NULL);
    //  memspace  = H5Screate_simple(4, chunk_dims, NULL);

    sprintf( str0, "/Z" );
    plist_id = H5Pcreate( H5P_DATASET_CREATE );
    H5Pset_chunk( plist_id, 4, chunk_dims );
    dset_id = H5Dcreate( file_id, str0, H5T_NATIVE_DOUBLE, filespace, H5P_DEFAULT, plist_id, H5P_DEFAULT );

    H5Pclose( plist_id );
    H5Sclose( filespace );

    filespace = H5Dget_space( dset_id );

    co = 0;

    PittPackReal Dz = ( F.Zb - F.Za ) / F.nChunk;

    for ( auto it = 0; it < F.nChunk; it++ )
    {
        // proc.enclosingBox(key,XYZ);
        Za = F.Za + it * Dz;
        Zb = F.Za + Dz * ( it + 1 );
        //      cout<<"Dz "<<Dz<<" Za "<<Za <<" Zb  "<<Zb<<endl;

        hz = N1 - 1.0;
        Zh = ( Zb - Za ) / ( hz );

        index = 0;

        for ( uint j = 0; j < L1; j++ )
        {
            for ( uint k = 0; k < M1; k++ )
            {
                for ( uint l = 0; l < N1; l++ )
                {
                    ztemp[index] = Za + Zh * l;
                    index++;
                }
            }
        }
        // define the offset, only in the fourth dimension

        offset[3] = myoffset + co;
        // printf("my_rank=%d offset=%d\n",my_rank,offset[3]);
        // status = H5Sselect_hyperslab(filespace, H5S_SELECT_SET, offset, NULL,
        // count, block);
        H5Sselect_hyperslab( filespace, H5S_SELECT_SET, offset, NULL, count, block );
        plist_id = H5Pcreate( H5P_DATASET_XFER );
        // H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);
        H5Pset_dxpl_mpio( plist_id, H5FD_MPIO_INDEPENDENT );
        // status = H5Dwrite(dset_id, H5T_NATIVE_DOUBLE,
        // memspace,filespace,plist_id, xtemp);
        if ( sizeof( PittPackReal ) == sizeof( double ) )
        {
            H5Dwrite( dset_id, H5T_NATIVE_DOUBLE, memspace, filespace, plist_id, ztemp );
        }
        else if ( sizeof( PittPackReal ) == sizeof( float ) )
        {
            H5Dwrite( dset_id, H5T_NATIVE_FLOAT, memspace, filespace, plist_id, ztemp );
        }
        H5Pclose( plist_id );
        co = co + 1;
    }
    //================================================================================
    //
    //                               write Q
    //
    //================================================================================
    // filespace = H5Screate_simple(4, dimsf, NULL);
    //  memspace  = H5Screate_simple(4, chunk_dims, NULL);
    dimsf[0] = F.nx;
    dimsf[1] = F.ny;
    dimsf[2] = F.nz;
    dimsf[3] = totalnumber;

    chunk_dims[0] = F.nx;
    chunk_dims[1] = F.ny;
    chunk_dims[2] = F.nz;
    chunk_dims[3] = 1;

    block[0] = chunk_dims[0];
    block[1] = chunk_dims[1];
    block[2] = chunk_dims[2];
    block[3] = chunk_dims[3];

    filespace = H5Screate_simple( 4, dimsf, NULL );
    memspace  = H5Screate_simple( 4, chunk_dims, NULL );

    sprintf( str0, "/Q" );
    plist_id = H5Pcreate( H5P_DATASET_CREATE );
    H5Pset_chunk( plist_id, 4, chunk_dims );
    dset_id = H5Dcreate( file_id, str0, H5T_NATIVE_DOUBLE, filespace, H5P_DEFAULT, plist_id, H5P_DEFAULT );

    H5Pclose( plist_id );
    H5Sclose( filespace );

    filespace = H5Dget_space( dset_id );

    // need to fix the offset for cell-center approach

    int my_rank, comsize;
    MPI_Comm_rank( MPI_COMM_WORLD, &my_rank );

    co = 0;

    for ( int it = 0; it < F.nChunk; it++ )
    {
        index = 0;

        for ( uint j = 0; j < L1 - 1; j++ )
        {
            for ( uint k = 0; k < M1 - 1; k++ )
            {
                for ( uint l = 0; l < N1 - 1; l++ )
                {
                    qtemp[index] = F( it, 0, j, k, l, 0 );
                    // cout<<" myrank "<<my_rank<< " qtemp "<<qtemp[index]<<endl;
                    index++;
                }
            }
        }
        // define the offset, only in the fourth dimension

        offset[3] = myoffset + co;
        // printf("my_rank=%d offset=%d\n",my_rank,offset[3]);
        // status = H5Sselect_hyperslab(filespace, H5S_SELECT_SET, offset, NULL,
        // count, block);
        H5Sselect_hyperslab( filespace, H5S_SELECT_SET, offset, NULL, count, block );
        plist_id = H5Pcreate( H5P_DATASET_XFER );
        // H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);
        H5Pset_dxpl_mpio( plist_id, H5FD_MPIO_INDEPENDENT );
        // status = H5Dwrite(dset_id, H5T_NATIVE_DOUBLE,
        // memspace,filespace,plist_id, xtemp);

        if ( sizeof( PittPackReal ) == sizeof( double ) )
        {
            H5Dwrite( dset_id, H5T_NATIVE_DOUBLE, memspace, filespace, plist_id, qtemp );
            //   cout<<BLUE "using double " RESET<<endl;
        }
        else if ( sizeof( PittPackReal ) == sizeof( float ) )
        {
            H5Dwrite( dset_id, H5T_NATIVE_FLOAT, memspace, filespace, plist_id, qtemp );
        }

        //      H5Dwrite(dset_id, H5T_NATIVE_DOUBLE, memspace,filespace,plist_id,
        // ztemp);
        H5Pclose( plist_id );
        co = co + 1;
    }

    // xdmfMultiBlock(F, com.mySize() ,com.myRank(), offset1,appx );
    xdmfMultiBlockCellCenter( F, com.mySize(), com.myRank(), offset1, appx );

    delete[] xtemp;
    delete[] ytemp;
    delete[] ztemp;
    delete[] qtemp;
}
#else

void Phdf5::writeMultiBlockCellCenter( ChunkedArray &F, uint appx, int dir, int aligndir )
{
    hid_t   file_id, dset_id;    /* file and dataset identifiers */
    hid_t   filespace, memspace; /* file and memory dataspace identifiers */
    hsize_t dimsf[4];            /* dataset dimensions */
    hsize_t chunk_dims[4];       /* chunk dimensions */
    hsize_t count[4];            /* hyperslab selection parameters */
    hsize_t block[4];
    hsize_t offset[4];
    hid_t   plist_id; /* property list identifier */
    uint    i, j, k, l;
    herr_t  status;
    // int         *data=NULL;

    /*
     * MPI variables
     */
    int           mpi_size, mpi_rank;
    PittPackReal *xtemp = NULL;
    PittPackReal *ytemp = NULL;
    PittPackReal *ztemp = NULL;
    PittPackReal *qtemp = NULL;

    uint L1 = F.nx + 1;
    uint M1 = F.ny + 1;
    uint N1 = F.nz + 1;

    xtemp = new PittPackReal[L1 * M1 * N1];
    ytemp = new PittPackReal[L1 * M1 * N1];
    ztemp = new PittPackReal[L1 * M1 * N1];
    qtemp = new PittPackReal[( L1 - 1 ) * ( M1 - 1 ) * ( N1 - 1 )];
    // cout<<L1<<"\t"<<M1<<"\t"<<N1<<endl;

    PittPackReal Xa;
    PittPackReal Xb;
    PittPackReal Ya;
    PittPackReal Yb;
    PittPackReal Za;
    PittPackReal Zb;

    char str[50];

    char str0[50];

    sprintf( str, H5FILE_NAME, appx );

    uint partialforestsize = F.nChunk; /*!<the forest size for each processor */

    CommPoint2Point<uint> com( &partialforestsize, 1 );
    uint                  offset1 = 0, totalvalue = 0;
    com.getOffset( partialforestsize, &offset1 );

    CommCollective<uint> comc( nullptr, 1 );
    comc.getTotalNumber( &offset1, &partialforestsize, &totalvalue );

    cout << " partial size " << partialforestsize << endl;

    hsize_t myoffset = offset1;
    //  total_size=totalvalue;

    totalnumber = totalvalue;
    /*
        cout << YELLOW "offsettttttttttttttt " << offset1 << RESET << endl;
        cout << YELLOW "L,M N " << L1 << " " << M1 << " " << N1 << RESET << endl;
        cout << YELLOW "totalvalue " << totalvalue << RESET << endl;
    */
    /*
     * Set up file access property list with parallel I/O access
     */

    plist_id = H5Pcreate( H5P_FILE_ACCESS );
    H5Pset_fapl_mpio( plist_id, MPI_COMM_WORLD, MPI_INFO_NULL );

    /*
     * Create a new file collectively and release property list identifier.
     */

    file_id = H5Fcreate( str, H5F_ACC_TRUNC, H5P_DEFAULT, plist_id );
    H5Pclose( plist_id );

    /*
     * Create the dataspace for the dataset.
     */
    /*
        dimsf[0] = L1;
        dimsf[1] = M1;
        dimsf[2] = N1;
        dimsf[3] = totalnumber;

        chunk_dims[0] = L1;
        chunk_dims[1] = M1;
        chunk_dims[2] = N1;
        chunk_dims[3] = 1;
    */
    dimsf[0] = N1;
    dimsf[1] = M1;
    dimsf[2] = L1;
    dimsf[3] = totalnumber;

    chunk_dims[0] = N1;
    chunk_dims[1] = M1;
    chunk_dims[2] = L1;
    chunk_dims[3] = 1;

    /*
     * Each process defines dataset in memory and writes it to the hyperslab
     * in the file.
     */

    count[0] = 1;
    count[1] = 1;
    count[2] = 1;
    count[3] = 1;

    block[0] = chunk_dims[0];
    block[1] = chunk_dims[1];
    block[2] = chunk_dims[2];
    block[3] = chunk_dims[3];

    /*
      for loop such that each processor can write all the blocks
    */
    offset[0] = 0;
    offset[1] = 0;
    offset[2] = 0;
    // offset[3] will be set at loop
    //================================================================================
    //
    //                               write X
    //
    //================================================================================
    // the first argument is the dimension which is 4

    filespace = H5Screate_simple( 4, dimsf, NULL );
    memspace  = H5Screate_simple( 4, chunk_dims, NULL );

    sprintf( str0, "/X" );
    plist_id = H5Pcreate( H5P_DATASET_CREATE );
    H5Pset_chunk( plist_id, 4, chunk_dims );
    dset_id = H5Dcreate( file_id, str0, H5T_NATIVE_DOUBLE, filespace, H5P_DEFAULT, plist_id, H5P_DEFAULT );

    H5Pclose( plist_id );
    H5Sclose( filespace );

    filespace = H5Dget_space( dset_id );

    unsigned int level, coord_index;

    PittPackReal dx, dy, dz;
    PittPackReal denum;
    PittPackReal XYZ[6];
    //  morton<N> key;
    hsize_t co = 0;

    PittPackReal Xh, Yh, Zh;
    PittPackReal hx, hy, hz;
    int          index;

    PittPackReal Dx = ( F.Xb - F.Xa ) / F.nChunk;

    PittPackReal Dy = ( F.Yb - F.Ya ) / F.nChunk;
    /*
        cout << YELLOW << "Xa " << F.Xa << " Xb " << F.Xb << RESET << endl;
        cout << YELLOW << "Ya " << F.Ya << " Yb " << F.Yb << RESET << endl;
        cout << YELLOW << "Za " << F.Za << " Zb " << F.Zb << RESET << endl;
        cout << YELLOW << "Xh " << ( F.Xb - F.Xa ) / double( L1 - 1. ) << RESET << endl;
        cout << YELLOW << "F.nChunk " << F.nChunk << " Dx " << Dx << RESET << endl;
    */
    int L[3] = {L1, M1, N1};

    for ( auto it = 0; it < F.nChunk; it++ )
    {
        // proc.enclosingBox(key,XYZ);

        if ( dir == 0 )
        {
            Xa = F.Xa + it * Dx;
            Xb = F.Xa + Dx * ( it + 1 );
            //  cout << BLUE << "Dx " << Dx << " Xa " << Xa << " Xb  " << Xb << RESET << endl;
        }
        else
        {
            Xa = F.Xa;
            Xb = F.Xb;
        }

        hx = L1 - 1.0;
        Xh = ( Xb - Xa ) / ( hx );
        /*
                index = 0;

                for ( uint j = 0; j < L1; j++ )
                {
                    for ( uint k = 0; k < M1; k++ )
                    {
                        for ( uint l = 0; l < N1; l++ )
                        {
                            xtemp[index] = Xa + Xh * j;
                            index++;
                        }
                    }
                }
        */

        getXcoord( L, Xa, Xh, aligndir, xtemp );

        // define the offset, only in the fourth dimension

        offset[3] = myoffset + co;

        printf( "offset=%llu\n", offset[3] );
        // status = H5Sselect_hyperslab(filespace, H5S_SELECT_SET, offset, NULL,
        // count, block);
        H5Sselect_hyperslab( filespace, H5S_SELECT_SET, offset, NULL, count, block );
        plist_id = H5Pcreate( H5P_DATASET_XFER );
        // H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);
        H5Pset_dxpl_mpio( plist_id, H5FD_MPIO_INDEPENDENT );
        // status = H5Dwrite(dset_id, H5T_NATIVE_DOUBLE,
        // memspace,filespace,plist_id, xtemp);

        if ( sizeof( PittPackReal ) == sizeof( double ) )
        {
            H5Dwrite( dset_id, H5T_NATIVE_DOUBLE, memspace, filespace, plist_id, xtemp );
        }
        else if ( sizeof( PittPackReal ) == sizeof( float ) )
        {
            H5Dwrite( dset_id, H5T_NATIVE_FLOAT, memspace, filespace, plist_id, xtemp );
        }

        H5Pclose( plist_id );
        co = co + 1;
    }

    H5Dclose( dset_id );

    // need to close dset_id since the new variable is going to be Y and Z ...
    // no need to open and close memspace and filespace all the time, simply reuse
    // them
    //================================================================================
    //
    //                               write Y
    //
    //================================================================================

    // filespace = H5Screate_simple(4, dimsf, NULL);
    //  memspace  = H5Screate_simple(4, chunk_dims, NULL);

    sprintf( str0, "/Y" );
    plist_id = H5Pcreate( H5P_DATASET_CREATE );
    H5Pset_chunk( plist_id, 4, chunk_dims );
    dset_id = H5Dcreate( file_id, str0, H5T_NATIVE_DOUBLE, filespace, H5P_DEFAULT, plist_id, H5P_DEFAULT );

    H5Pclose( plist_id );
    H5Sclose( filespace );

    filespace = H5Dget_space( dset_id );

    co = 0;

    for ( auto it = 0; it < F.nChunk; it++ )
    {
        //       proc.enclosingBox(key,XYZ);
        //       cout<<XYZ[0]<<" "<<XYZ[1]<<" "<<XYZ[2]<<" "<<XYZ[3]<<" "<<XYZ[4]<<"
        // "<<XYZ[5]<<endl;
        if ( dir == 1 )
        {
            Ya = F.Ya + it * Dy;
            Yb = F.Ya + Dy * ( it + 1 );
        }
        else
        {
            Ya = F.Ya;
            Yb = F.Yb;
        }

        hy = M1 - 1.0;
        Yh = ( Yb - Ya ) / ( hy );

        index = 0;
        /*
                for ( uint j = 0; j < L1; j++ )
                {
                    for ( uint k = 0; k < M1; k++ )
                    {
                        for ( uint l = 0; l < N1; l++ )
                        {
                            ytemp[index] = Ya + Yh * k;
                            index++;
                        }
                    }
                }
                // define the offset, only in the fourth dimension
        */
        getYcoord( L, Ya, Yh, aligndir, ytemp );

        offset[3] = myoffset + co;
        // printf("my_rank=%d offset=%d\n",my_rank,offset[3]);
        // status = H5Sselect_hyperslab(filespace, H5S_SELECT_SET, offset, NULL,
        // count, block);
        H5Sselect_hyperslab( filespace, H5S_SELECT_SET, offset, NULL, count, block );
        plist_id = H5Pcreate( H5P_DATASET_XFER );
        // H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);
        // status = H5Dwrite(dset_id, H5T_NATIVE_DOUBLE,
        // memspace,filespace,plist_id, xtemp);
        if ( sizeof( PittPackReal ) == sizeof( double ) )
        {
            H5Dwrite( dset_id, H5T_NATIVE_DOUBLE, memspace, filespace, plist_id, ytemp );
        }
        else if ( sizeof( PittPackReal ) == sizeof( float ) )
        {
            H5Dwrite( dset_id, H5T_NATIVE_FLOAT, memspace, filespace, plist_id, ytemp );
        }

        H5Pclose( plist_id );
        co = co + 1;
    }

    //================================================================================
    //
    //                               write Z
    //
    //================================================================================

    // filespace = H5Screate_simple(4, dimsf, NULL);
    //  memspace  = H5Screate_simple(4, chunk_dims, NULL);

    sprintf( str0, "/Z" );
    plist_id = H5Pcreate( H5P_DATASET_CREATE );
    H5Pset_chunk( plist_id, 4, chunk_dims );
    dset_id = H5Dcreate( file_id, str0, H5T_NATIVE_DOUBLE, filespace, H5P_DEFAULT, plist_id, H5P_DEFAULT );

    H5Pclose( plist_id );
    H5Sclose( filespace );

    filespace = H5Dget_space( dset_id );

    co = 0;

    PittPackReal Dz = ( F.Zb - F.Za ) / F.nChunk;

    for ( auto it = 0; it < F.nChunk; it++ )
    {
        if ( dir == 2 )
        {
            Za = F.Za + it * Dz;
            Zb = F.Za + Dz * ( it + 1 );
        }
        else
        {
            Za = F.Za;
            Zb = F.Zb;
        }
        hz = N1 - 1.0;
        Zh = ( Zb - Za ) / ( hz );
        /*
                index = 0;

                for ( uint j = 0; j < L1; j++ )
                {
                    for ( uint k = 0; k < M1; k++ )
                    {
                        for ( uint l = 0; l < N1; l++ )
                        {
                            ztemp[index] = Za + Zh * l;
                            index++;
                        }
                    }
                }

        */
        // define the offset, only in the fourth dimension

        getZcoord( L, Za, Zh, aligndir, ztemp );

        offset[3] = myoffset + co;
        // printf("my_rank=%d offset=%d\n",my_rank,offset[3]);
        // status = H5Sselect_hyperslab(filespace, H5S_SELECT_SET, offset, NULL,
        // count, block);
        H5Sselect_hyperslab( filespace, H5S_SELECT_SET, offset, NULL, count, block );
        plist_id = H5Pcreate( H5P_DATASET_XFER );
        // H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);
        H5Pset_dxpl_mpio( plist_id, H5FD_MPIO_INDEPENDENT );
        // status = H5Dwrite(dset_id, H5T_NATIVE_DOUBLE,
        // memspace,filespace,plist_id, xtemp);
        if ( sizeof( PittPackReal ) == sizeof( double ) )
        {
            H5Dwrite( dset_id, H5T_NATIVE_DOUBLE, memspace, filespace, plist_id, ztemp );
        }
        else if ( sizeof( PittPackReal ) == sizeof( float ) )
        {
            H5Dwrite( dset_id, H5T_NATIVE_FLOAT, memspace, filespace, plist_id, ztemp );
        }
        H5Pclose( plist_id );
        co = co + 1;
    }
    //================================================================================
    //
    //                               write Q
    //
    //================================================================================
    // filespace = H5Screate_simple(4, dimsf, NULL);
    //  memspace  = H5Screate_simple(4, chunk_dims, NULL);
    // my z corresponds to hdf5 x and vice versa
    /*
        dimsf[0] = F.nx;
        dimsf[1] = F.ny;
        dimsf[2] = F.nz;
        dimsf[3] = totalnumber;

        chunk_dims[0] = F.nx;
        chunk_dims[1] = F.ny;
        chunk_dims[2] = F.nz;
        chunk_dims[3] = 1;
    */

    dimsf[0] = F.nz;
    dimsf[1] = F.ny;
    dimsf[2] = F.nx;
    dimsf[3] = totalnumber;

    chunk_dims[0] = F.nz;
    chunk_dims[1] = F.ny;
    chunk_dims[2] = F.nx;
    chunk_dims[3] = 1;

    block[0] = chunk_dims[0];
    block[1] = chunk_dims[1];
    block[2] = chunk_dims[2];
    block[3] = chunk_dims[3];

    filespace = H5Screate_simple( 4, dimsf, NULL );
    memspace  = H5Screate_simple( 4, chunk_dims, NULL );

    sprintf( str0, "/Q" );
    plist_id = H5Pcreate( H5P_DATASET_CREATE );
    H5Pset_chunk( plist_id, 4, chunk_dims );
    dset_id = H5Dcreate( file_id, str0, H5T_NATIVE_DOUBLE, filespace, H5P_DEFAULT, plist_id, H5P_DEFAULT );

    H5Pclose( plist_id );
    H5Sclose( filespace );

    filespace = H5Dget_space( dset_id );

    // need to fix the offset for cell-center approach

    int my_rank, comsize;
    MPI_Comm_rank( MPI_COMM_WORLD, &my_rank );

    co = 0;

    for ( int it = 0; it < F.nChunk; it++ )
    {
        /*
                index = 0;

                for ( uint j = 0; j < L1 - 1; j++ )
                {
                    for ( uint k = 0; k < M1 - 1; k++ )
                    {
                        for ( uint l = 0; l < N1 - 1; l++ )
                        {
                            qtemp[index] = F( it, dir, j, k, l, 0 );
                            //                cout<<" myrank "<<my_rank<< " qtemp "<<qtemp[index]<<endl;
                            index++;
                        }
                    }
                }
        */

        getQ( F, it, L, aligndir, qtemp );
        // define the offset, only in the fourth dimension

        offset[3] = myoffset + co;
        // printf("my_rank=%d offset=%d\n",my_rank,offset[3]);
        // status = H5Sselect_hyperslab(filespace, H5S_SELECT_SET, offset, NULL,
        // count, block);
        H5Sselect_hyperslab( filespace, H5S_SELECT_SET, offset, NULL, count, block );
        plist_id = H5Pcreate( H5P_DATASET_XFER );
        // H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);
        H5Pset_dxpl_mpio( plist_id, H5FD_MPIO_INDEPENDENT );
        // status = H5Dwrite(dset_id, H5T_NATIVE_DOUBLE,
        // memspace,filespace,plist_id, xtemp);

        if ( sizeof( PittPackReal ) == sizeof( double ) )
        {
            H5Dwrite( dset_id, H5T_NATIVE_DOUBLE, memspace, filespace, plist_id, qtemp );
            //   cout<<BLUE "using double " RESET<<endl;
        }
        else if ( sizeof( PittPackReal ) == sizeof( float ) )
        {
            H5Dwrite( dset_id, H5T_NATIVE_FLOAT, memspace, filespace, plist_id, qtemp );
        }

        //      H5Dwrite(dset_id, H5T_NATIVE_DOUBLE, memspace,filespace,plist_id,
        // ztemp);
        H5Pclose( plist_id );
        co = co + 1;
    }

    // xdmfMultiBlock(F, com.mySize() ,com.myRank(), offset1,appx );
    xdmfMultiBlockCellCenter( F, com.mySize(), com.myRank(), offset1, appx );

    delete[] xtemp;
    delete[] ytemp;
    delete[] ztemp;
    delete[] qtemp;
}

#endif

void Phdf5::xdmfMultiBlockCellCenter( ChunkedArray &F, integer comsize, integer my_rank, uint offset, uint appx )
{
    MPI_File     fp;
    int          buf[1000], np = comsize;
    MPI_Request  request;
    unsigned int i;
    int          j;

    MPI_Status  status;
    const char *names[] = {"X", "Y", "Z"};

    // chunks in hdf5 work in reverse order my x corresponds to hdf5 z, and vice verse
    uint L1 = F.nz + 1;
    uint M1 = F.ny + 1;
    uint N1 = F.nx + 1;

    char strname[1000], fname[1000];
    uint ncube_total = totalnumber;
    sprintf( strname, XDMF_NAME, appx );
    char str[1000];

    int offset0 = 156;
    int offset1 = 2005;
    //  const int offset1 = 1992;

    cout << "size str" << strlen( strname ) << endl;

    MPI_File_open( MPI_COMM_WORLD, strname, MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &fp );

    // printf("%s\n",str);
    // dont need to specify the address here as they are both reside in one
    // folder, that is
    // why I use H5FILE in addition to H5FILE_NAME
    sprintf( fname, H5FILE, appx );
    char strLCell[100];
    char strMCell[100];
    char strNCell[100];

    char strL[100];
    char strM[100];
    char strN[100];
    char strNcube[1000];
    char stroff[1000];
    int  index;

    if ( L1 < 10 )
    {
        sprintf( strL, "00%d", L1 );
        sprintf( strM, "00%d", M1 );
        sprintf( strN, "00%d", N1 );

        sprintf( strLCell, "00%d", L1 - 1 );
        sprintf( strMCell, "00%d", M1 - 1 );
        sprintf( strNCell, "00%d", N1 - 1 );
    }
    else if ( L1 < 100 )
    {
        sprintf( strL, "0%d", L1 );
        sprintf( strM, "0%d", M1 );
        sprintf( strN, "0%d", N1 );

        sprintf( strLCell, "0%d", L1 - 1 );
        sprintf( strMCell, "0%d", M1 - 1 );
        sprintf( strNCell, "0%d", N1 - 1 );
    }
    else if ( L1 < 1000 )
    {
        sprintf( strL, "%d", L1 );
        sprintf( strM, "%d", M1 );
        sprintf( strN, "%d", N1 );

        sprintf( strLCell, "%d", L1 - 1 );
        sprintf( strMCell, "%d", M1 - 1 );
        sprintf( strNCell, "%d", N1 - 1 );
    }

    else
    {
        printf( "discretization too big go change xdmf function\n" );
        exit( 0 );
    }

    if ( ncube_total < 10 )
    {
        sprintf( strNcube, "00000%d", ncube_total );
    }
    else if ( ncube_total < 100 )
    {
        sprintf( strNcube, "0000%d", ncube_total );
    }
    else if ( ncube_total < 1000 )
    {
        sprintf( strNcube, "000%d", ncube_total );
    }
    else if ( ncube_total < 10000 )
    {
        sprintf( strNcube, "00%d", ncube_total );
    }
    else if ( ncube_total < 100000 )
    {
        sprintf( strNcube, "0%d", ncube_total );
    }
    else if ( ncube_total < 1000000 )
    {
        sprintf( strNcube, "%d", ncube_total );
    }
    else
    {
        printf( "number of cubes larger than 1000000\n" );
        exit( 0 );
    }

    int  ab[2];
    uint co = 0;
    // MPI_Status status;
    // MPI_Request request;

    if ( my_rank != 0 )
    {
        MPI_Irecv( ab, 2, MPI_INT, 0, 0, MPI_COMM_WORLD, &request );
    }

    // a counts the offset for header whic is only written by process rank 0 and
    // and b the hyperslab part for each cube
    int a = 0, b = 0;
    if ( my_rank == 0 )
    {
        sprintf( str, "<?xml version=\"1.0\" ?>\n" );
        MPI_File_write( fp, str, strlen( str ), MPI_CHAR, &status );
        a = a + strlen( str );
        sprintf( str, "<!DOCTYPE Xdmf SYSTEM \"Xdmf.dtd\" []> \n" );
        MPI_File_write( fp, str, strlen( str ), MPI_CHAR, &status );
        a = a + strlen( str );
        sprintf( str, "<Xdmf Version=\"2.0\">\n" );
        MPI_File_write( fp, str, strlen( str ), MPI_CHAR, &status );
        a = a + strlen( str );
        sprintf( str, "<Domain>\n" );
        MPI_File_write( fp, str, strlen( str ), MPI_CHAR, &status );
        a = a + strlen( str );
        sprintf( str, "<Grid Name=\"AMR\" GridType=\"Collection\" "
                      "CollectionType=\"Spatial\">\n" );
        MPI_File_write( fp, str, strlen( str ), MPI_CHAR, &status );
        a = a + strlen( str );
        printf( "size a=%d\n", a );
        // for(unsigned int i=0;i<cube.size();i++)
        co = 0;

        // for(auto it=F.trees.begin();it!=F.trees.end();it++)
        {
            //        for(auto it2=(*it).begin();it2!=(*it).end();it2++)
            for ( auto it2 = 0; it2 < F.nChunk; it2++ )
            {
                // b calculates the offset for each block, need to be set to zero
                // becasue it is inside the loop
                // but we dont want to accumulate, calculate only for one block
                b = 0;
                sprintf( str, "   <Grid Name=\"mesh0\" GridType=\"Uniform\">\n" );
                MPI_File_write( fp, str, strlen( str ), MPI_CHAR, &status );
                b = b + strlen( str );
                sprintf( str,
                         "        <Topology TopologyType=\"3DSMesh\" "
                         "NumberOfElements=\"%s %s %s\"/>\n",
                         strL, strM, strN );
                // printf("hyperslab =%d\n",strlen(str));
                MPI_File_write( fp, str, strlen( str ), MPI_CHAR, &status );
                b = b + strlen( str );
                sprintf( str, "          <Geometry GeometryType=\"X_Y_Z\">  \n" );
                MPI_File_write( fp, str, strlen( str ), MPI_CHAR, &status );
                b = b + strlen( str );

                for ( j = 0; j < 3; j++ )
                {
                    sprintf( str,
                             "      <DataItem ItemType=\"HyperSlab\" Dimensions=\"%s "
                             "%s %s %d\" NumberType=\"Float\" Precision=\"4\"  "
                             "Type=\"HyperSlab\"> \n",
                             strL, strM, strN, 1 );
                    MPI_File_write( fp, str, strlen( str ), MPI_CHAR, &status );
                    b = b + strlen( str );
                    sprintf( str, "         <DataItem Dimensions=\"3 4\" Format=\"XML\" > \n" );
                    MPI_File_write( fp, str, strlen( str ), MPI_CHAR, &status );
                    b     = b + strlen( str );
                    index = offset + co;
                    integer_string( stroff, index );
                    sprintf( str, "         %d %d %d %s  \n", 0, 0, 0, stroff );
                    MPI_File_write( fp, str, strlen( str ), MPI_CHAR, &status );
                    b = b + strlen( str );
                    sprintf( str, "         %d %d %d %d  \n", 1, 1, 1, 1 );
                    MPI_File_write( fp, str, strlen( str ), MPI_CHAR, &status );
                    b = b + strlen( str );
                    sprintf( str, "         %s %s %s %d  \n", strL, strM, strN, 1 );
                    MPI_File_write( fp, str, strlen( str ), MPI_CHAR, &status );
                    b = b + strlen( str );
                    sprintf( str, "         </DataItem>\n" );
                    MPI_File_write( fp, str, strlen( str ), MPI_CHAR, &status );
                    b = b + strlen( str );
                    sprintf( str,
                             "          <DataItem Name=\"%s\" Dimensions=\"%s %s %s "
                             "%s\" NumberType=\"Float\" Precision=\"4\" "
                             "Format=\"HDF\">\n",
                             names[j], strL, strM, strN, strNcube );
                    MPI_File_write( fp, str, strlen( str ), MPI_CHAR, &status );
                    b = b + strlen( str );
                    sprintf( str, "          %s:/%s\n", fname, names[j] );
                    MPI_File_write( fp, str, strlen( str ), MPI_CHAR, &status );
                    b = b + strlen( str );
                    sprintf( str, "         </DataItem>\n" );
                    MPI_File_write( fp, str, strlen( str ), MPI_CHAR, &status );
                    b = b + strlen( str );
                    sprintf( str, "         </DataItem>\n" );
                    MPI_File_write( fp, str, strlen( str ), MPI_CHAR, &status );
                    b = b + strlen( str );
                }
                sprintf( str, "      </Geometry>   \n" );
                MPI_File_write( fp, str, strlen( str ), MPI_CHAR, &status );
                b = b + strlen( str );
                // added because of including solution vector Q
                sprintf( str, "         <Attribute Name=\"Q\" AttributeType=\"Scalar\" "
                              "Center=\"Cell\"> \n" );
                MPI_File_write( fp, str, strlen( str ), MPI_CHAR, &status );
                b = b + strlen( str );
                sprintf( str,
                         "      <DataItem ItemType=\"HyperSlab\" Dimensions=\"%s "
                         "%s %s %d\" NumberType=\"Float\" Precision=\"4\"  "
                         "Type=\"HyperSlab\"> \n",
                         strLCell, strMCell, strNCell, 1 );
                MPI_File_write( fp, str, strlen( str ), MPI_CHAR, &status );
                b = b + strlen( str );
                sprintf( str, "         <DataItem Dimensions=\"3 4\" Format=\"XML\" > \n" );
                MPI_File_write( fp, str, strlen( str ), MPI_CHAR, &status );
                b     = b + strlen( str );
                index = offset + co;
                integer_string( stroff, index );
                sprintf( str, "         %d %d %d %s  \n", 0, 0, 0, stroff );
                MPI_File_write( fp, str, strlen( str ), MPI_CHAR, &status );
                b = b + strlen( str );
                sprintf( str, "         %d %d %d %d  \n", 1, 1, 1, 1 );
                MPI_File_write( fp, str, strlen( str ), MPI_CHAR, &status );
                b = b + strlen( str );
                sprintf( str, "         %s %s %s %d  \n", strLCell, strMCell, strNCell, 1 );
                MPI_File_write( fp, str, strlen( str ), MPI_CHAR, &status );
                b = b + strlen( str );
                sprintf( str, "         </DataItem>\n" );
                MPI_File_write( fp, str, strlen( str ), MPI_CHAR, &status );
                b = b + strlen( str );
                sprintf( str,
                         "          <DataItem Name=\"Q\" Dimensions=\"%s %s %s "
                         "%s\" NumberType=\"Float\" Precision=\"4\" "
                         "Format=\"HDF\">\n",
                         strLCell, strMCell, strNCell, strNcube );
                MPI_File_write( fp, str, strlen( str ), MPI_CHAR, &status );
                b = b + strlen( str );
                sprintf( str, "          %s:/Q\n", fname );
                MPI_File_write( fp, str, strlen( str ), MPI_CHAR, &status );
                b = b + strlen( str );

                sprintf( str, "         </DataItem>\n" );
                MPI_File_write( fp, str, strlen( str ), MPI_CHAR, &status );
                b = b + strlen( str );
                sprintf( str, "         </DataItem>\n" );
                MPI_File_write( fp, str, strlen( str ), MPI_CHAR, &status );
                b = b + strlen( str );
                sprintf( str, " </Attribute>\n" );
                MPI_File_write( fp, str, strlen( str ), MPI_CHAR, &status );
                b = b + strlen( str );
                sprintf( str, "  </Grid>\n" );
                MPI_File_write( fp, str, strlen( str ), MPI_CHAR, &status );
                b = b + strlen( str );
#if ( 1 )
                co++;
/*
                // put error here
                if ( a != offset0 || b != offset1 )
                {
                    printf( "??????????????????????????????????????????\n" );
                    printf( "Go fix your offset for Xdmf meta data file a=%d offset0=%d "
                            "b=%d offset1=%d\n",
                            a, offset0, b, offset1 );
                    printf( "??????????????????????????????????????????\n" );
//                    exit( 0 );
                }
*/
//              printf("Go fix your offset for Xdmf meta data file a=%d
// offset0=%d b=%d offset1=%d\n",a,offset0,b,offset1);
#endif
            }
        }

        ab[0] = a;
        ab[1] = b;
        // send the values of a and b to other procs
        for ( int i = 1; i < np; i++ )
        {
            MPI_Send( ab, 2, MPI_INT, i, 0, MPI_COMM_WORLD );
        }
    }
#if ( 1 )
    else
    {
        //   calcOff(&offset0,offset1);
        // other ranks do something else
        // 5 is added because of the header of the file written by  processor 0
        MPI_Wait( &request, &status );
        offset0 = ab[0];
        offset1 = ab[1];

        // printf("offset0(%d) = %d  offset1(%d) = %d\n",my_rank,offset0,my_rank,offset1);
        int mpi_offset = 0;
        co             = 0;

        for ( int it = 0; it < F.nChunk; it++ )
        {
            mpi_offset = ( ( offset + co ) * offset1 + offset0 );
            mpi_offset = mpi_offset * sizeof( char );
            // printf("my_rank = %d offset =%d mpi_offset=%d
            // i=%d\n",my_rank,offset,mpi_offset,i);

            MPI_File_seek( fp, mpi_offset, MPI_SEEK_SET );
            // MPI_File_set_view(fp,mpi_offset,MPI_CHAR,MPI_CHAR, "native",
            // MPI_INFO_NULL );
            sprintf( str, "   <Grid Name=\"mesh0\" GridType=\"Uniform\">\n" );
            MPI_File_write( fp, str, strlen( str ), MPI_CHAR, &status );
            sprintf( str,
                     "        <Topology TopologyType=\"3DSMesh\" "
                     "NumberOfElements=\"%s %s %s\"/>\n",
                     strL, strM, strN );
            MPI_File_write( fp, str, strlen( str ), MPI_CHAR, &status );
            sprintf( str, "          <Geometry GeometryType=\"X_Y_Z\">  \n" );
            MPI_File_write( fp, str, strlen( str ), MPI_CHAR, &status );

#if ( 1 )
            for ( j = 0; j < 3; j++ )
            {
                sprintf( str,
                         "      <DataItem ItemType=\"HyperSlab\" Dimensions=\"%s "
                         "%s %s %d\" NumberType=\"Float\" Precision=\"4\"  "
                         "Type=\"HyperSlab\"> \n",
                         strL, strM, strN, 1 );
                MPI_File_write( fp, str, strlen( str ), MPI_CHAR, &status );
                sprintf( str, "         <DataItem Dimensions=\"3 4\" Format=\"XML\" > \n" );
                MPI_File_write( fp, str, strlen( str ), MPI_CHAR, &status );
                index = offset + co;
                integer_string( stroff, index );
                sprintf( str, "         %d %d %d %s  \n", 0, 0, 0, stroff );
                MPI_File_write( fp, str, strlen( str ), MPI_CHAR, &status );
                sprintf( str, "         %d %d %d %d  \n", 1, 1, 1, 1 );
                MPI_File_write( fp, str, strlen( str ), MPI_CHAR, &status );
                sprintf( str, "         %s %s %s %d  \n", strL, strM, strN, 1 );
                MPI_File_write( fp, str, strlen( str ), MPI_CHAR, &status );
                sprintf( str, "         </DataItem>\n" );
                MPI_File_write( fp, str, strlen( str ), MPI_CHAR, &status );
                sprintf( str,
                         "          <DataItem Name=\"%s\" Dimensions=\"%s %s %s "
                         "%s\" NumberType=\"Float\" Precision=\"4\" "
                         "Format=\"HDF\">\n",
                         names[j], strL, strM, strN, strNcube );
                MPI_File_write( fp, str, strlen( str ), MPI_CHAR, &status );
                sprintf( str, "          %s:/%s\n", fname, names[j] );
                MPI_File_write( fp, str, strlen( str ), MPI_CHAR, &status );
                sprintf( str, "         </DataItem>\n" );
                MPI_File_write( fp, str, strlen( str ), MPI_CHAR, &status );
                sprintf( str, "         </DataItem>\n" );
                MPI_File_write( fp, str, strlen( str ), MPI_CHAR, &status );
            }
#endif
            sprintf( str, "      </Geometry>   \n" );
            MPI_File_write( fp, str, strlen( str ), MPI_CHAR, &status );

// insert here the value of qs
#if ( 1 )
            sprintf( str, "         <Attribute Name=\"Q\" AttributeType=\"Scalar\" "
                          "Center=\"Cell\"> \n" );
            MPI_File_write( fp, str, strlen( str ), MPI_CHAR, &status );
            sprintf( str,
                     "      <DataItem ItemType=\"HyperSlab\" Dimensions=\"%s %s "
                     "%s %d\" NumberType=\"Float\" Precision=\"4\"  "
                     "Type=\"HyperSlab\"> \n",
                     strLCell, strMCell, strNCell, 1 );
            MPI_File_write( fp, str, strlen( str ), MPI_CHAR, &status );
            sprintf( str, "         <DataItem Dimensions=\"3 4\" Format=\"XML\" > \n" );
            MPI_File_write( fp, str, strlen( str ), MPI_CHAR, &status );
            index = offset + co;
            integer_string( stroff, index );
            sprintf( str, "         %d %d %d %s  \n", 0, 0, 0, stroff );
            MPI_File_write( fp, str, strlen( str ), MPI_CHAR, &status );
            sprintf( str, "         %d %d %d %d  \n", 1, 1, 1, 1 );
            MPI_File_write( fp, str, strlen( str ), MPI_CHAR, &status );
            sprintf( str, "         %s %s %s %d  \n", strLCell, strMCell, strNCell, 1 );
            MPI_File_write( fp, str, strlen( str ), MPI_CHAR, &status );
            sprintf( str, "         </DataItem>\n" );
            MPI_File_write( fp, str, strlen( str ), MPI_CHAR, &status );
            sprintf( str,
                     "          <DataItem Name=\"Q\" Dimensions=\"%s %s %s %s\" "
                     "NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n",
                     strLCell, strMCell, strNCell, strNcube );
            MPI_File_write( fp, str, strlen( str ), MPI_CHAR, &status );
            sprintf( str, "          %s:/Q\n", fname );
            MPI_File_write( fp, str, strlen( str ), MPI_CHAR, &status );
            sprintf( str, "         </DataItem>\n" );
            MPI_File_write( fp, str, strlen( str ), MPI_CHAR, &status );
            sprintf( str, "         </DataItem>\n" );
            MPI_File_write( fp, str, strlen( str ), MPI_CHAR, &status );
            sprintf( str, " </Attribute>\n" );
            MPI_File_write( fp, str, strlen( str ), MPI_CHAR, &status );

            sprintf( str, "  </Grid>\n" );
            MPI_File_write( fp, str, strlen( str ), MPI_CHAR, &status );
            co++;
#endif
        }
        // add the last part of the xdmf text file to close the arguments
    }
#endif
    MPI_Barrier( MPI_COMM_WORLD );
#if ( 1 )
    if ( my_rank == ( np - 1 ) )
    {
        sprintf( str, "  </Grid>\n" );
        MPI_File_write( fp, str, strlen( str ), MPI_CHAR, &status );
        sprintf( str, "      </Domain>   \n" );
        MPI_File_write( fp, str, strlen( str ), MPI_CHAR, &status );
        sprintf( str, "  </Xdmf>\n" );
        MPI_File_write( fp, str, strlen( str ), MPI_CHAR, &status );
    }
#endif
    MPI_File_sync( fp );
    MPI_File_close( &fp );
}

void Phdf5::getXcoord( int *L, const double Xa, const double Xh, const int aligndir, double *xtemp )
{
    int index = 0;
    int L1    = L[0];
    int M1    = L[1];
    int N1    = L[2];
    // this is consistent with the major directions specfied in rearrange
    // n0(0,0,1) slices are in xy-plane
    // major direction=1 is planar rotation in y irection but still slices are in xy plane
    // n2(1,0,0)

    // cout<< " ============================= "<< " aligndir "<< aligndir  <<endl;
    if ( aligndir == 0 )
    {
        for ( uint k = 0; k < N1; k++ )
        {
            for ( uint j = 0; j < M1; j++ )
            {
                for ( uint i = 0; i < L1; i++ )
                {
                    xtemp[index] = Xa + Xh * i;
                    // cout<< " xtemp "<<xtemp[index]<<endl;
                    index++;
                }
            }
        }
    }
    else if ( aligndir == 1 )
    {
        for ( uint k = 0; k < N1; k++ )
        {
            for ( uint i = 0; i < L1; i++ )
            {
                for ( uint j = 0; j < M1; j++ )
                {
                    {
                        xtemp[index] = Xa + Xh * i;
                        index++;
                    }
                }
            }
        }
    }
    else if ( aligndir == 2 )
    {
        for ( uint i = 0; i < N1; i++ )

        {
            for ( uint j = 0; j < M1; j++ )
            {
                for ( uint k = 0; k < L1; k++ )
                {
                    xtemp[index] = Xa + Xh * i;
                    index++;
                }
            }
        }
    }
    else
    {
        cout << " undefined direction in phdf5" << endl;
        exit( 0 );
    }
    //  cout<< " ============================= "<<endl;
}

void Phdf5::getYcoord( int *L, const double Ya, const double Yh, const int aligndir, double *ytemp )
{
    int index = 0;
    int L1    = L[0];
    int M1    = L[1];
    int N1    = L[2];

    if ( aligndir == 0 )
    {
        for ( uint k = 0; k < N1; k++ )
        {
            for ( uint j = 0; j < M1; j++ )
            {
                for ( uint i = 0; i < L1; i++ )
                {
                    ytemp[index] = Ya + Yh * j;
                    //  cout<< " ytemp "<<ytemp[index]<<endl;
                    index++;
                }
            }
        }
    }
    else if ( aligndir == 1 )
    {
        for ( uint k = 0; k < N1; k++ )
        {
            for ( uint i = 0; i < L1; i++ )
            {
                for ( uint j = 0; j < M1; j++ )

                {
                    ytemp[index] = Ya + Yh * j;
                    index++;
                }
            }
        }
    }
    else if ( aligndir == 2 )
    {
        for ( uint i = 0; i < L1; i++ )
        {
            for ( uint j = 0; j < M1; j++ )

            {
                for ( uint k = 0; k < N1; k++ )
                {
                    ytemp[index] = Ya + Yh * j;
                    index++;
                }
            }
        }
    }
    else
    {
        cout << " undefined direction in phdf5" << endl;
        exit( 0 );
    }
    //  cout<< " ============================= "<<endl;
}

void Phdf5::getZcoord( int *L, const double Za, const double Zh, const int aligndir, double *ztemp )
{
    int index = 0;
    int L1    = L[0];
    int M1    = L[1];
    int N1    = L[2];

    // z is the same for two directions 0 and 1 since rottion in on n(0,0,0) plane

    if ( aligndir == 0 )
    {
        for ( uint k = 0; k < N1; k++ )
        {
            for ( uint j = 0; j < M1; j++ )
            {
                for ( uint i = 0; i < L1; i++ )
                {
                    ztemp[index] = Za + Zh * k;
                    //  cout<< " ztemp "<<ztemp[index]<<endl;
                    index++;
                }
            }
        }
    }
    else if ( aligndir == 1 )
    {
        for ( uint k = 0; k < N1; k++ )
        {
            for ( uint i = 0; i < L1; i++ )
            {
                for ( uint j = 0; j < M1; j++ )

                {
                    ztemp[index] = Za + Zh * k;
                    index++;
                }
            }
        }
    }
    else if ( aligndir == 2 )
    {
        for ( uint i = 0; i < L1; i++ )
        {
            for ( uint j = 0; j < M1; j++ )

            {
                for ( uint k = 0; k < N1; k++ )
                {
                    ztemp[index] = Za + Zh * k;
                    index++;
                }
            }
        }
    }
    else
    {
        cout << " undefined direction in phdf5" << endl;
        exit( 0 );
    }

    cout << " ============================= " << endl;
}

void Phdf5::getQ( ChunkedArray &F, const int chunkId, int *L, const int aligndir, double *qtemp )
{
    int index = 0;
    int L1    = L[0] - 1;
    int M1    = L[1] - 1;
    int N1    = L[2] - 1;
    // this is consistent with the major directions specfied in rearrange
    // n0(0,0,1) slices are in xy-plane
    // major direction=1 is planar rotation in y irection but still slices are in xy plane
    // n2(1,0,0)

    //    cout << "Q ============================= "
    //         << " aligndir " << aligndir << endl;
    if ( aligndir == 0 )
    {
        for ( uint k = 0; k < N1; k++ )
        {
            for ( uint j = 0; j < M1; j++ )
            {
                for ( uint i = 0; i < L1; i++ )
                {
                    qtemp[index] = F( chunkId, aligndir, i, j, k, 0 );
                    //    cout<< " qtemp "<<qtemp[index]<<endl;
                    index++;
                }
            }
        }
    }

    else if ( aligndir == 1 )
    {
        for ( uint k = 0; k < N1; k++ )
        {
            for ( uint i = 0; i < L1; i++ )
            {
                for ( uint j = 0; j < M1; j++ )
                {
                    {
                        qtemp[index] = F( chunkId, aligndir, i, j, k, 0 );
                        index++;
                    }
                }
            }
        }
    }
    else if ( aligndir == 2 )
    {
        for ( uint i = 0; i < N1; i++ )

        {
            for ( uint j = 0; j < M1; j++ )
            {
                for ( uint k = 0; k < L1; k++ )
                {
                    qtemp[index] = F( chunkId, aligndir, i, j, k, 0 );
                    index++;
                }
            }
        }
    }
    else
    {
        cout << " undefined direction in phdf5" << endl;
        exit( 0 );
    }

    // cout << " ============================= " << endl;
}
