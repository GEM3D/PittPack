#include "triDiag.h"
#include "mathFunction.h"
#include "params.h"

void TriDiag::setElems( int nCh, int nzCh, double *sub, double *sup )
{
    nChunk  = nCh;
    nzChunk = nzCh;

    subDiag = new double[3];
    supDiag = new double[3];

    for ( int i = 0; i < 3; i++ )
    {
        subDiag[i] = sub[i];
        supDiag[i] = sup[i];
    }

#pragma acc enter data create( this [0:1] )
#pragma acc update device( this )
#pragma acc enter data create( subDiag [0:3] )
#pragma acc update device( subDiag [0:3] )
#pragma acc enter data create( supDiag [0:3] )
#pragma acc update device( supDiag [0:3] )
}

TriDiag::~TriDiag()
{
    if ( subDiag != NULL )
    {
        delete[] subDiag;
    }
    if ( supDiag != NULL )
    {
        delete[] supDiag;
    }

#if ( OPENACC )
#pragma acc exit data delete ( subDiag )
#pragma acc exit data delete ( supDiag )
#pragma acc exit data delete ( this )
#endif
}

#if ( OPENACC )
#pragma acc routine
#endif
void TriDiag::thomas( ChunkedArray &P, double *onDiag, const int i, const int j, const int dir, const int index )
{
    double bet;

    int result;
    // the enterior is always set according to eigenvalues
    // the two ends decided by BC

    double Sn[ZSIZE];
    int    this_rank = 1;

#if ( 1 )
    bet                         = onDiag[0];
    P( 0, dir, i, j, 0, index ) = P( 0, dir, i, j, 0, index ) / ( bet );

    if ( absolute( bet ) < 1e-12 )
    {
        result = THOMAS_FAIL;
    }

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
    bet   = onDiag[1] - subDiag[1] * Sn[k];

    if ( absolute( bet ) < 1e-12 )
    {
        result = THOMAS_FAIL;
    }
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
    int kend   = nzChunk;
    int count  = 2;

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
            bet       = onDiag[1] - subDiag[1] * Sn[count];

            if ( absolute( bet ) < 1e-12 )
            {
                result = THOMAS_FAIL;
            }

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

    k     = nChunk * nzChunk - 1;
    Sn[k] = supDiag[1] / bet;
    bet   = onDiag[2] - subDiag[2] * Sn[k];
    if ( absolute( bet ) < 1e-12 )
    {
        result = THOMAS_FAIL;
    }

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
        // if ( id != 0 ) this is a bug, fixed now, backsun goes all the way back
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
    // return(result);
}

#if ( OPENACC )
#pragma acc routine
#endif
void TriDiag::thomasSingleBlock( ChunkedArray &P, double *onDiag, int i, int j, int dir, int index )
{
    double bet;

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
    //  int Nx = nxChunk;
    //  int Ny = nyChunk;

    //    double onDiag[3];
    //    onDiag[1] = getEigenVal( i, j );

#if ( 1 )
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
    Sn[k]                       = supDiag[0] / bet;
    bet                         = onDiag[1] - subDiag[1] * Sn[k];
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
        Sn[k]                       = supDiag[1] / bet;
        bet                         = onDiag[1] - subDiag[1] * Sn[k];
        P( 0, dir, i, j, k, index ) = ( P( 0, dir, i, j, k, index ) - subDiag[1] * P( 0, dir, i, j, k - 1, index ) ) / bet;
    }

    k = nzChunk - 1;

    Sn[k]                       = supDiag[1] / bet;
    bet                         = onDiag[2] - subDiag[2] * Sn[k];
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

#if ( 1 )
#if ( OPENACC )
#pragma acc routine
#endif
void TriDiag::thomasPeriodic( ChunkedArray &P, double *onDiag, int i, int j, int dir, int index )
{
    double bet;

    // the enterior is always set according to eigenvalues
    // the two ends decided by BC

    //    double onDiag[3];

    double alpha = 1.0;
    double beta  = 1.0;
    double Sn[ZSIZE];
    double Zn[ZSIZE];

    // onDiag[1] = -2.0 + ( -2.0 + 2. * cosine( ( i + num[0] ) * pi / ( Nx + denum[0] ) ) )
    //          + ( -2.0 + 2. * cosine( ( j + num[1] ) * pi / ( Ny + denum[1] ) ) );

    //    onDiag[1] = getEigenVal( i, j );

    int this_rank = 0;
    int myRank    = 0;

    // assign bc
    // bc 0 >> dirichlet, value known at ghost point

    onDiag[0] = onDiag[1];
    onDiag[2] = onDiag[1];

    double gamma = -onDiag[0];

    onDiag[0] = onDiag[0] - gamma;
    onDiag[2] = onDiag[2] - alpha * beta / gamma;

#if ( 1 )
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
    Sn[k]                       = supDiag[0] / bet;
    bet                         = onDiag[1] - subDiag[1] * Sn[k];
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
    int kend   = nzChunk;
    int count  = 2;

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
            Sn[count]                    = supDiag[1] / bet;
            bet                          = onDiag[1] - subDiag[1] * Sn[count];
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

    k     = nChunk * nzChunk - 1;
    Sn[k] = supDiag[1] / bet;
    bet   = onDiag[2] - subDiag[2] * Sn[k];
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

#endif

    count = 0;
    for ( int id = 0; id < nChunk; id++ )
    {
        for ( int k = 0; k < nzChunk; k++ )
        {
            Zn[count]                    = P( id, dir, i, j, k, index );
            P( id, dir, i, j, k, index ) = 0.0;
            count++;
        }
    }

    P( 0, 0, i, j, 0, index )                    = gamma;
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
    Sn[k]                       = supDiag[0] / bet;
    bet                         = onDiag[1] - subDiag[1] * Sn[k];
    P( 0, dir, i, j, k, index ) = ( P( 0, dir, i, j, k, index ) - subDiag[1] * P( 0, dir, i, j, k - 1, index ) ) / bet;

#if ( DEBUG )
    if ( myRank == this_rank )
    {
        cout << RED << "my_rank " << myRank << " Pn [ " << 1 << "] " << P( 0, dir, i, j, k, index ) << RESET << endl;
    }
#endif
    kstart = 2;
    kend   = nzChunk;
    count  = 2;

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
            Sn[count]                    = supDiag[1] / bet;
            bet                          = onDiag[1] - subDiag[1] * Sn[count];
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

    k     = nChunk * nzChunk - 1;
    Sn[k] = supDiag[1] / bet;
    bet   = onDiag[2] - subDiag[2] * Sn[k];
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
}
#endif

#if ( OPENACC )
#pragma acc routine
#endif
int TriDiag::thomasReal( ChunkedArray &P, double *onDiag, const int i, const int j, const int dir )
{
    int    index = 0;
    double bet;
    int    this_rank = 1;
    int    result    = SUCCESS;

#if ( DEBUG )
    for ( int id = 0; id < nChunk; id++ )
    {
        for ( int k = 0; k < nzChunk; k++ )
        {
            // cout<<RED<<"my_rank "<<myRank<< " nz = "<<nzChunk<<" "<<P(id,0,i,j,k,0)<<RESET<<endl;
            cout << RED << "Thomas before for rank "
                 << " " << P( id, dir, i, j, k, index ) << RESET << endl;
        }
    }
#endif

#if ( 1 )
    bet = onDiag[0];

    if ( absolute( bet ) < 1e-12 )
    {
        result = THOMAS_FAIL;
    }

    P( 0, dir, i, j, 0, index ) = P( 0, dir, i, j, 0, index ) / ( bet );

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

    // Sn[k] = supDiag[0] / bet;

    P( 0, dir, i, j, k, 1 ) = supDiag[0] / bet;
    bet                     = onDiag[1] - subDiag[1] * P( 0, dir, i, j, k, 1 );

    if ( absolute( bet ) < 1e-12 )
    {
        result = THOMAS_FAIL;
    }

    P( 0, dir, i, j, k, index ) = ( P( 0, dir, i, j, k, index ) - subDiag[1] * P( 0, dir, i, j, k - 1, index ) ) / bet;

#if ( DEBUG )
    if ( myRank == this_rank )
    {
        cout << RED << "my_rank " << myRank << " Pn [ " << 1 << "] " << P( 0, dir, i, j, k, index ) << RESET << endl;
    }
#endif
    int kstart = 2;
    int kend   = nzChunk;
    int count  = 2;

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
            // Sn[count] = supDiag[1] / bet;
            P( id, dir, i, j, k, 1 ) = supDiag[1] / bet;
            // bet = onDiag[1] - subDiag[1] * Sn[count];
            bet = onDiag[1] - subDiag[1] * P( id, dir, i, j, k, 1 );

            if ( absolute( bet ) < 1e-12 )
            {
                result = THOMAS_FAIL;
            }

            P( id, dir, i, j, k, index ) = ( P( id, dir, i, j, k, index ) - subDiag[1] * P( id, dir, i, j, k - 1, index ) ) / bet;
#if ( DEBUG )
            if ( myRank == this_rank )
            {
                cout << RED << "my_rank " << myRank << " Pn [ " << count << "] " << P( id, dir, i, j, k, index ) << RESET << endl;
            }
#endif
            //            count++;
        }
    }

    //   k = nChunk * nzChunk - 1;
    // Sn[k] = supDiag[1] / bet;
    P( nChunk - 1, dir, i, j, ( nzChunk - 1 ), 1 ) = supDiag[1] / bet;

    //    bet = onDiag[2] - subDiag[2] * Sn[k];
    bet = onDiag[2] - subDiag[2] * P( nChunk - 1, dir, i, j, ( nzChunk - 1 ), 1 );

    if ( absolute( bet ) < 1e-12 )
    {
        result = THOMAS_FAIL;
    }

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
#if ( 1 )
    count = nChunk * nzChunk - 2;

    for ( int id = nChunk - 1; id >= 0; id-- )
    {
        //  if ( id != 0 )
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
            // P( id, dir, i, j, k, index ) -= Sn[count + 1] * P( id, dir, i, j, k + 1, index );
            P( id, dir, i, j, k, index ) -= P( id, dir, i, j, k + 1, 1 ) * P( id, dir, i, j, k + 1, index );
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
            cout << RED << "Thomas result for rank "
                 << " " << P( id, dir, i, j, k, index ) << RESET << endl;
        }
    }
#endif

    /*
    for(int id=0;id<nChunk;id++)
    {
    for(int k=0;k<nzChunk;k++)
    {
    //   P(id,0,i,j,k,0)=1000.0;
       P(id,0,i,j,k,1)=0.0;
    }
    }
    */

    // return(THOMAS_FAIL);
    return ( result );
}

#if ( OPENACC )
#pragma acc routine
#endif
void TriDiag::thomasPeriodicReal( ChunkedArray &P, double *onDiag, int i, int j, int dir )
{
    double bet;
    int    index = 0;

    // the enterior is always set according to eigenvalues
    // the two ends decided by BC
    //    double onDiag[3];

    double alpha = 1.0;
    double beta  = 1.0;
    double Zn[ZSIZE];

    // onDiag[1] = -2.0 + ( -2.0 + 2. * cosine( ( i + num[0] ) * pi / ( Nx + denum[0] ) ) )
    //          + ( -2.0 + 2. * cosine( ( j + num[1] ) * pi / ( Ny + denum[1] ) ) );

    //    onDiag[1] = getEigenVal( i, j );

    int this_rank = 0;
    int myRank    = 0;

    // assign bc
    // bc 0 >> dirichlet, value known at ghost point

    onDiag[0] = onDiag[1];
    onDiag[2] = onDiag[1];

    double gamma = -onDiag[0];

    onDiag[0] = onDiag[0] - gamma;
    onDiag[2] = onDiag[2] - alpha * beta / gamma;

#if ( 1 )
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
    //  Sn[k] = supDiag[0] / bet;
    //  bet = onDiag[1] - subDiag[1] * Sn[k];

    P( 0, dir, i, j, k, 1 )     = supDiag[0] / bet;
    bet                         = onDiag[1] - subDiag[1] * P( 0, dir, i, j, k, 1 );
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
    int kend   = nzChunk;
    int count  = 2;

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
            P( id, dir, i, j, k, 1 ) = supDiag[1] / bet;
            bet                      = onDiag[1] - subDiag[1] * P( id, dir, i, j, k, 1 );

            P( id, dir, i, j, k, index ) = ( P( id, dir, i, j, k, index ) - subDiag[1] * P( id, dir, i, j, k - 1, index ) ) / bet;
            //   count++;
        }
    }

    P( nChunk - 1, dir, i, j, ( nzChunk - 1 ), 1 ) = supDiag[1] / bet;
    bet                                            = onDiag[2] - subDiag[2] * P( nChunk - 1, dir, i, j, ( nzChunk - 1 ), 1 );
    P( nChunk - 1, dir, i, j, ( nzChunk - 1 ), index )
    = ( P( nChunk - 1, dir, i, j, ( nzChunk - 1 ), index ) - subDiag[2] * P( nChunk - 1, dir, i, j, nzChunk - 2, index ) ) / bet;

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
            //    P( id, dir, i, j, k, index ) -= Sn[count + 1] * P( id, dir, i, j, k + 1, index );
            P( id, dir, i, j, k, index ) -= P( id, dir, i, j, k + 1, 1 ) * P( id, dir, i, j, k + 1, index );

            if ( myRank == this_rank )
            {
                //                cout << BLUE << "my_rank " << myRank << " Pn[   " << count << " ]" << P( id, 0, i, j, k, 0 ) << RESET <<
                // endl;
            }
            count--;
        }
    }
#endif

#endif

    count = 0;
    for ( int id = 0; id < nChunk; id++ )
    {
        for ( int k = 0; k < nzChunk; k++ )
        {
            Zn[count]                    = P( id, dir, i, j, k, 0 );
            P( id, dir, i, j, k, index ) = 0.0;
            count++;
        }
    }

    P( 0, 0, i, j, 0, index )                    = gamma;
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
    P( 0, dir, i, j, k, 1 ) = supDiag[0] / bet;
    bet                     = onDiag[1] - subDiag[1] * P( 0, dir, i, j, k, 1 );

    //    Sn[k] = supDiag[0] / bet;
    //    bet = onDiag[1] - subDiag[1] * Sn[k];
    P( 0, dir, i, j, k, index ) = ( P( 0, dir, i, j, k, index ) - subDiag[1] * P( 0, dir, i, j, k - 1, index ) ) / bet;

#if ( DEBUG )
    if ( myRank == this_rank )
    {
        cout << RED << "my_rank " << myRank << " Pn [ " << 1 << "] " << P( 0, dir, i, j, k, index ) << RESET << endl;
    }
#endif
    kstart = 2;
    kend   = nzChunk;
    count  = 2;

    for ( int id = 0; id < nChunk; id++ )
    {
        // if ( id != 0 ) this was a bug as the solution has to calculate till the first element
        {
            kstart = 0;
        }
        if ( id == ( nChunk - 1 ) )
        {
            kend = nzChunk - 1;
        }

        for ( int k = kstart; k < kend; k++ )
        {
            P( id, dir, i, j, k, 1 )     = supDiag[1] / bet;
            bet                          = onDiag[1] - subDiag[1] * P( id, dir, i, j, k, 1 );
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

    P( nChunk - 1, dir, i, j, ( nzChunk - 1 ), 1 ) = supDiag[1] / bet;
    bet                                            = onDiag[2] - subDiag[2] * P( nChunk - 1, dir, i, j, ( nzChunk - 1 ), 1 );
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
            P( id, dir, i, j, k, index ) -= P( id, dir, i, j, k + 1, 1 ) * P( id, dir, i, j, k + 1, index );

            // P( id, dir, i, j, k, index ) -= Sn[count + 1] * P( id, dir, i, j, k + 1, index );

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
}
