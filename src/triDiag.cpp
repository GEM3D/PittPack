#include "triDiag.hpp"
#include "mathFunction.hpp"
#include "params.h"

void TriDiag::setElems( int nCh, int nzCh, double *sub, double *sup )
{
    nChunk  = nCh;
    nzChunk = nzCh;

    subDiag = new double[3];
    supDiag = new double[3];
    bc      = new char[2];

    /*
        bc[0]=BC[4];
        bc[1]=BC[5];

    */

    for ( int i = 0; i < 3; i++ )
    {
        subDiag[i] = sub[i];
        supDiag[i] = sup[i];
    }

#if ( PITTPACKACC )
#pragma acc enter data create( this [0:1] )
#pragma acc update device( this )
#pragma acc enter data create( subDiag [0:3] )
#pragma acc update device( subDiag [0:3] )
#pragma acc enter data create( supDiag [0:3] )
#pragma acc update device( supDiag [0:3] )
#pragma acc enter data create( bc [0:2] )
#endif
}

void TriDiag::assignBC( char *BC )
{
    bc[0] = BC[4];
    bc[1] = BC[5];

#if ( PITTPACKACC )
#pragma acc update device( bc [0:2] )
#endif
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
    if ( bc != NULL )
    {
        delete[] bc;
    }

#if ( PITTPACKACC )
#pragma acc exit data delete ( subDiag )
#pragma acc exit data delete ( supDiag )
#pragma acc exit data delete ( bc )
#pragma acc exit data delete ( this )
#endif
}

#if ( PITTPACKACC )
#pragma acc routine
#endif
void TriDiag::thomas( ChunkedArray &P, double *onDiag, const int i, const int j, const int dir, const int index )
{
    double bet;

    int result;
    // the enterior is always set according to eigenvalues
    // the two ends decided by BC

    double Sn[ZSIZE];
    // int    this_rank = 1;

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

#if ( PITTPACKACC )
#pragma acc routine
#endif
void TriDiag::thomasSingleBlock( ChunkedArray &P, double *onDiag, int i, int j, int dir, int index )
{
    double bet;

    double Sn[ZSIZE];
    // the enterior is always set according to eigenvalues
    // the two ends decided by BC
    /*
    #if(PITTPACKACC)
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
#if ( PITTPACKACC )
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

#if ( PITTPACKACC )
#pragma acc routine
#endif
int TriDiag::thomasReal( ChunkedArray &P, double *onDiag, const int i, const int j, const int dir )
{
    int    index = 0;
    double bet;
    // int    this_rank = 1;
    int result = SUCCESS;

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

#if ( PITTPACKACC )
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

// only works for 2^n-1
#if ( PITTPACKACC )
#pragma acc routine vector
#endif
void TriDiag::crp( const int n, double offdiag, double *tmpA, double *tmpC, double *tmpRHS, double *thmA, double *thmC, double *thmB,
                   double *gam1, double *rhs )
{
    int level = myLog2( n );
    int d0;
    int d1;
    int m = n;
    // int d02;

    double a0, c0, cte; // aMinus1, cMinus1, aPlus1, cPlus1, cte;
    int    ind0;
    int    indMinus1;
    int    indPlus1;

    tmpA[0] = ( offdiag );
    tmpC[0] = ( offdiag );

#if ( PITTPACKACC )
#pragma acc loop seq
#endif
    for ( int i = 0; i < level - 1; i++ )
    {
        d1 = 1 << ( i + 1 );
        d0 = 1 << i;
        m  = m >> 1;

        a0 = tmpA[i];
        c0 = tmpC[i];

        cte = 1. - ( 2. * a0 * c0 );
        //  printf("%lf %lf\n",a0,c0);
#if ( 1 )
#if ( PITTPACKACC )
#pragma acc loop vector firstprivate( cte ) private( ind0, indMinus1, indPlus1 )
#endif
        //#pragma acc loop vector
        for ( int j = 0; j < m; j = j + 1 )
        {
            ind0      = d1 * j + ( d1 - 1 );
            indMinus1 = d1 * j + ( d1 - 1 ) - d1 / 2;
            indPlus1  = d1 * j + ( d1 - 1 ) + d1 / 2;

            rhs[d1 * j + ( d1 - 1 )] = ( rhs[indMinus1] * ( -a0 ) + rhs[ind0] + rhs[indPlus1] * ( -c0 ) ) / cte;
            // rhs[d1 * j + (d1 - 1)] = (rhs[indMinus1]* (-a0) +  rhs[ind0] + rhs[indPlus1]*(-c0));

            //     printf("index %d uses %d %d %d a0=%lf c0=%lf cte=%lf rhs=%lf \n ",j,ind0,indMinus1,indPlus1,a0,c0,cte,rhs[d1*j+(d1-1)]);
        }
#endif
        // printf("============================\n");
        tmpA[i + 1] = -a0 * a0 / cte;
        tmpC[i + 1] = -c0 * c0 / cte;
    }

    /*
      printf("before back sub \n");
      printf("rhs=\n");
      for (int i = 0; i < n; i++) {
        printf("%lf\n", rhs[i]);
      }

      printf("tmpA = \n");
      for (int i = 0; i < level+1; i++) {
        printf("tmpA = %lf tmpC=%lf\n", tmpA[i],tmpC[i]);
      }
    */
    //
    // solution of a three by three matrix

    // printf("m =%d \n",m);

    m = 3;
#if ( PITTPACKACC )
#pragma acc loop seq
#endif
    for ( int i = 0; i < 1; i++ )
    {
        d1        = 1 << ( level );
        ind0      = ( d1 - 1 );
        indMinus1 = ( d1 - 1 ) - d1 / 2;
        indPlus1  = ( d1 - 1 ) + d1 / 2;

        tmpRHS[1] = rhs[ind0];
        tmpRHS[0] = rhs[indMinus1];
        tmpRHS[2] = rhs[indPlus1];

        thmA[1] = tmpA[level - 1];
        thmA[2] = tmpA[level - 1];
        thmC[0] = tmpC[level - 1];
        thmC[1] = tmpC[level - 1];
        //    double thmB[3]={1.0,1.0,1.0};
        //    double gam1[3];

#if ( 1 )
        thomasLowMem( 3, thmA, thmB, thmC, tmpRHS, gam1 );

        rhs[ind0]      = tmpRHS[1];
        rhs[indMinus1] = tmpRHS[0];
        rhs[indPlus1]  = tmpRHS[2];
    }
    //   printf("tmpRHS= %lf %lf %lf\n",tmpRHS[0],tmpRHS[1],tmpRHS[2]);
    // calculate the backsub stage
    int div;

    // breaking data dependency so that we can perform the tasks in parallel
#if ( PITTPACKACC )
#pragma acc loop seq
#endif
    for ( int i = level - 1; i > 0; i-- )
    {
        // d0=-2.0;
        d1 = 1 << ( i - 1 );
        d0 = 1 << i;
//     printf("   i =%d d1 =  %d m =%d d0=%d\n ",i, d1,m,d0);

//  printf("--------------------m= %d div =%d ---------------\n",m,div);

//#pragma acc loop vector
#if ( PITTPACKACC )
#pragma acc loop vector firstprivate( d0, m )
#endif
        //#pragma acc loop seq
        for ( int j = 0; j < m; j = j + 1 )
        {
            rhs[d0 * j + ( d0 - 1 ) - d0 / 2] = rhs[d0 * j + ( d0 - 1 ) - d0 / 2] - tmpA[i - 1] * rhs[d0 * j + ( d0 - 1 )];
            //     printf("tmpA[]=%lf\n",tmpA[i]);
            //    printf("index(%d) uses %d  previous %d rhs = %lf b= %lf middle %lf\n ",j,
            //    d0*j+(d0-1),d0*j+(d0-1)-d0/2,rhs[d0*j+(d0-1)-d0/2],b[d0*j+(d0-1)-d0/2],b[d0*j+(d0-1)]  );
        }

//#pragma acc loop vector
#if ( PITTPACKACC )
#pragma acc loop vector firstprivate( d0, m )
#endif
        //#pragma acc loop seq
        for ( int j = 0; j < m; j = j + 1 )
        {
            rhs[d0 * j + ( d0 - 1 ) + d0 / 2] = rhs[d0 * j + ( d0 - 1 ) + d0 / 2] - tmpC[i - 1] * rhs[d0 * j + ( d0 - 1 )];
            //   printf("index(%d) uses %d  previous %d rhs = %lf b= %lf middle %lf \n ",j,
            //   d0*j+(d0-1),d0*j+(d0-1)+d0/2,rhs[d0*j+(d0-1)+d0/2], b[d0*j+(d0-1)+d0/2],b[d0*j+(d0-1)]  );
        }

        // m = 2 * ((level) - i + 1) + 1;

        div = 1 << ( level - i + 1 ) + 1;

        //    m = div/2 + 1;
        m = div - 1;

        //  printf("--------------------m= %d div =%d ---------------\n",m,div);
    }

#endif
}

#if ( PITTPACKACC )
#pragma acc routine vector
#endif
void TriDiag::pcr( int n, double *a, double *c, double *d )
{
    int level = myLog2( n ) + 1;
    int s;

    double r;
    // these go on register
    int m = VECLENGTH;
    // int size = n / m + 1;

    double a0[20];
    double a1[20];
    double a2[20];

    double c0[20];
    double c1[20];
    double c2[20];

    double d0[20];
    double d1[20];
    double d2[20];

    int index;

#if ( PITTPACKACC )
#pragma acc loop seq private( a0, a1, a2, c0, c1, c2, d0, d1, d2 )
#endif
    for ( int p = 0; p < level; p++ )
    {
        s = 1 << p;

#if ( PITTPACKACC )
#pragma acc loop private( a0, a1, a2, c0, c1, c2, d0, d1, d2, index )
#endif
        for ( int i = 0; i < n; i++ )
        {
            index = i / m;

            if ( i - s < 0 && i + s < n )
            {
                a0[index] = a[i];
                a1[index] = 0.0;
                a2[index] = a[i + s];

                c0[index] = c[i];
                c1[index] = 0.0;
                c2[index] = c[i + s];

                d0[index] = d[i];
                d1[index] = 0.0;
                d2[index] = d[i + s];
            }
            else if ( i + s >= n && i - s >= 0 )
            {
                a0[index] = a[i];
                a1[index] = a[i - s];
                a2[index] = 0.0;

                c0[index] = c[i];
                c1[index] = c[i - s];
                c2[index] = 0.0;

                d0[index] = d[i];
                d1[index] = d[i - s];
                d2[index] = 0.0;
            }
            // both indices are ok to assign
            else
            {
                a0[index] = a[i];
                a1[index] = a[i - s];
                a2[index] = a[i + s];

                c0[index] = c[i];
                c1[index] = c[i - s];
                c2[index] = c[i + s];

                d0[index] = d[i];
                d1[index] = d[i - s];
                d2[index] = d[i + s];
            }
        }

        // both indices are ok to assign
        //

#if ( PITTPACKACC )
#pragma acc loop private( a0, a1, a2, c0, c1, c2, d0, d1, d2, r, index )
#endif
        for ( int i = 0; i < n; i++ )
        {
            index = i / m;
            r     = 1. / ( 1. - a0[index] * c1[index] - c0[index] * a2[index] );
            a[i]  = -r * a0[index] * a1[index];
            c[i]  = -r * c0[index] * c2[index];
            d[i]  = r * ( d0[index] - a0[index] * d1[index] - c0[index] * d2[index] );
        }

        //    printf("i+s = %d  i+s =  %d\n",i-s,i+s );
    }
}

#if ( PITTPACKACC )
#pragma acc routine seq
#endif
bool TriDiag::checkRhs( const double *rh, const int N )
{
    bool   bl  = 1;
    double tol = 1.e-10;

    for ( int i = 0; i < N; i++ )
    {
        if ( rh[i] > tol )
        {
            bl = 0;
            break;
        }
    }
    return ( bl );
}
#if ( PITTPACKACC )
#pragma acc routine seq
#endif
void TriDiag::thomasLowMemNoBCV2( double *tmpMG, double *rh, double *diag, int index )
{
    double bet;

   

    //   double a[3],c[3];
    double b[3];

    int N = nChunk * nzChunk;
    /*
     #if(PITTPACKACC)
     #pragma acc loop seq
     #endif
         for ( int i = 0; i < 3; i++ )
         {
             a[i] = subDiag[i];
             c[i] = supDiag[i];
         }
     */
    b[0] = diag[0];
    b[1] = diag[1];
    b[2] = diag[2];

    rh[0] = rh[0] / ( bet = b[0] );

    int j    = 1;
    //tmpMG[j] = supDiag[j - 1] / bet;
    //periodic only, supDiag[0]=0.0
    tmpMG[j] = 0.0;
    bet      = b[1] - subDiag[1] * tmpMG[j];
    rh[1]    = ( rh[1] - subDiag[1] * rh[j - 1] ) / bet;

#if ( PITTPACKACC )
#pragma acc loop seq
#endif
    for ( int j = 2; j < N - 1; j++ )
    {
        //       DecompositioN and forward substitution.
        tmpMG[j] = supDiag[1] / bet;
        bet      = b[1] - subDiag[1] * tmpMG[j];
        rh[j]    = ( rh[j] - subDiag[1] * rh[j - 1] ) / bet;
    }

    j        = N - 1;
    tmpMG[j] = supDiag[1] / bet;
    
    bet      = b[2] -0.0* subDiag[2] * tmpMG[j];
    rh[j]    = ( rh[j] - 0.*subDiag[2] * rh[j - 1] ) / bet;

    //  cout << a[2] << " " << b[2] << eNdl;
    //  cout<<RED<<rh[0]<<RESET<<endl;
#if ( PITTPACKACC )
#pragma acc loop seq
#endif
    for ( int j = ( N - 2 ); j >= 0; j-- )
    {
        rh[j] -= tmpMG[j + 1] * rh[j + 1];
    }

    // cout<<RED<<rh[0]<<RESET<<endl;
}

#if ( PITTPACKACC )
#pragma acc routine seq
#endif
void TriDiag::thomasLowMemNoBCV1( double *tmpMG, double *rh, double *diag, int index )
{
    double bet;

   

    //   double a[3],c[3];
    double b[3];

    int N = nChunk * nzChunk;
    /*
     #if(PITTPACKACC)
     #pragma acc loop seq
     #endif
         for ( int i = 0; i < 3; i++ )
         {
             a[i] = subDiag[i];
             c[i] = supDiag[i];
         }
     */
    b[0] = diag[0];
    b[1] = diag[1];
    b[2] = diag[2];

    rh[0] = rh[0] / ( bet = b[0] );

    int j    = 1;
    //tmpMG[j] = supDiag[j - 1] / bet;
    //periodic only, supDiag[0]=0.0
    tmpMG[j] = 0.0;
    bet      = b[1] - subDiag[1] * tmpMG[j];
    rh[1]    = ( rh[1] - subDiag[1] * rh[j - 1] ) / bet;

#if ( PITTPACKACC )
#pragma acc loop seq
#endif
    for ( int j = 2; j < N - 1; j++ )
    {
        //       DecompositioN and forward substitution.
        tmpMG[j] = supDiag[1] / bet;
        bet      = b[1] - subDiag[1] * tmpMG[j];
        rh[j]    = ( rh[j] - subDiag[1] * rh[j - 1] ) / bet;
    }

    j        = N - 1;
    tmpMG[j] = supDiag[1] / bet;
    bet      = b[2] - subDiag[2] * tmpMG[j];
    rh[j]    = ( rh[j] - subDiag[2] * rh[j - 1] ) / bet;

    //  cout << a[2] << " " << b[2] << eNdl;
    //  cout<<RED<<rh[0]<<RESET<<endl;
#if ( PITTPACKACC )
#pragma acc loop seq
#endif
    for ( int j = ( N - 2 ); j >= 0; j-- )
    {
        rh[j] -= tmpMG[j + 1] * rh[j + 1];
    }

    // cout<<RED<<rh[0]<<RESET<<endl;
}
//  no boundary version

#if ( PITTPACKACC )
#pragma acc routine seq
#endif
void TriDiag::thomasLowMem( double *tmpMG, double *rh, double diag, int index )
{
    double bet;

    //   double a[3],c[3];
    double b[3];

    int n = nChunk * nzChunk;
    int N = n;
    /*
    #if(PITTPACKACC)
    #pragma acc loop seq
    #endif
        for ( int i = 0; i < 3; i++ )
        {
            a[i] = subDiag[i];
            c[i] = supDiag[i];
        }
    */
    b[0] = diag;
    b[1] = diag;
    b[2] = diag;

    // this inserted to prevent NN-NN-NN from blowing up
    if ( fabs( diag + 2. ) < 1.e-10 )
    {
        return;
    }

    // for Dirirchlet
    if ( bc[0] == 'D' )
    {
        b[0] = b[1] - 1.;
    }
    if ( bc[1] == 'D' )
    {
        b[2] = b[1] - 1.;
    }
    // for Neumann, all modifications are done on the stencil

    if ( bc[0] == 'N' )
    {
        b[0] = b[1] + 1.;
    }
    if ( bc[1] == 'N' )
    {
        b[2] = b[1] + 1.0;
    }

    rh[0] = rh[0] / ( bet = b[0] );

    int j    = 1;
    tmpMG[j] = supDiag[j - 1] / bet;
    bet      = b[1] - subDiag[1] * tmpMG[j];
    rh[1]    = ( rh[1] - subDiag[1] * rh[j - 1] ) / bet;

#if ( PITTPACKACC )
#pragma acc loop seq
#endif
    for ( int j = 2; j < n - 1; j++ )
    {
        //       DecompositioN and forward substitution.
        tmpMG[j] = supDiag[1] / bet;
        bet      = b[1] - subDiag[1] * tmpMG[j];
        rh[j]    = ( rh[j] - subDiag[1] * rh[j - 1] ) / bet;
    }

    j        = N - 1;
    tmpMG[j] = supDiag[1] / bet;
    bet      = b[2] - subDiag[2] * tmpMG[j];
    rh[j]    = ( rh[j] - subDiag[2] * rh[j - 1] ) / bet;

    //  cout << a[2] << " " << b[2] << eNdl;
#if ( PITTPACKACC )
#pragma acc loop seq
#endif
    for ( int j = ( n - 2 ); j >= 0; j-- )
    {
        rh[j] -= tmpMG[j + 1] * rh[j + 1];
        // cout << " j " << j << eNdl;
    }
}


//  no boundary version
#if ( PITTPACKACC )
#pragma acc routine seq
#endif
void TriDiag::thomasLowMemNoBC( double *tmpMG, double *rh, double *diag, int index )
{
    double bet;

   

    //   double a[3],c[3];
    double b[3];

    int N = nChunk * nzChunk;
    /*
     #if(PITTPACKACC)
     #pragma acc loop seq
     #endif
         for ( int i = 0; i < 3; i++ )
         {
             a[i] = subDiag[i];
             c[i] = supDiag[i];
         }
     */
    b[0] = diag[0];
    b[1] = diag[1];
    b[2] = diag[2];

    rh[0] = rh[0] / ( bet = b[0] );

    int j    = 1;
    tmpMG[j] = supDiag[j - 1] / bet;
    //periodic only, supDiag[0]=0.0
    //tmpMG[j] = 0.0;
    bet      = b[1] - subDiag[1] * tmpMG[j];
    rh[1]    = ( rh[1] - subDiag[1] * rh[j - 1] ) / bet;

#if ( PITTPACKACC )
#pragma acc loop seq
#endif
    for ( int j = 2; j < N - 1; j++ )
    {
        //       DecompositioN and forward substitution.
        tmpMG[j] = supDiag[1] / bet;
        bet      = b[1] - subDiag[1] * tmpMG[j];
        rh[j]    = ( rh[j] - subDiag[1] * rh[j - 1] ) / bet;
    }

    j        = N - 1;
    tmpMG[j] = supDiag[1] / bet;
    bet      = b[2] - subDiag[2] * tmpMG[j];
    rh[j]    = ( rh[j] - subDiag[2] * rh[j - 1] ) / bet;

    //  cout << a[2] << " " << b[2] << eNdl;
    //  cout<<RED<<rh[0]<<RESET<<endl;
#if ( PITTPACKACC )
#pragma acc loop seq
#endif
    for ( int j = ( N - 2 ); j >= 0; j-- )
    {
        rh[j] -= tmpMG[j + 1] * rh[j + 1];
    }

    // cout<<RED<<rh[0]<<RESET<<endl;
}
//  no boundary version
#if ( PITTPACKACC )
#pragma acc routine seq
#endif
void TriDiag::thomasLowMem( int N, double *a, double *b, double *c, double *r, double *gam )
{
    double bet;
    int    n = N;

    /*
        if ( b[1] == 0.0 )
        {
            cout<<b[0]<<endl;
            throw std::runtime_error( "Error in Tridag" );
        }
    */
    // If this happens then you shoul/d rewrite your equations as a set of order N
    // 1 ,with u 2 trivially eliminated.

    r[0] = r[0] / ( bet = b[0] );
    // cout << r[0] << endl;

    int j  = 1;
    gam[j] = c[j - 1] / bet;
    bet    = b[1] - a[1] * gam[j];
    r[1]   = ( r[1] - a[1] * r[j - 1] ) / bet;

#if ( PITTPACKACC )
#pragma acc loop seq
#endif
    for ( int j = 2; j < n - 1; j++ )
    {
        //       Decomposition and forward substitution.
        gam[j] = c[1] / bet;
        bet    = b[1] - a[1] * gam[j];
        r[j]   = ( r[j] - a[1] * r[j - 1] ) / bet;
    }

    j      = N - 1;
    gam[j] = c[1] / bet;
    bet    = b[2] - a[2] * gam[j];
    r[j]   = ( r[j] - a[2] * r[j - 1] ) / bet;

    // cout << a[2] << " " << b[2] << endl;
/*
  for (int j = 0; j < n; j++)
  {
    cout << "u_low_mem " << r[j] << endl;
  }
*/
#if ( PITTPACKACC )
#pragma acc loop seq
#endif
    for ( int j = ( n - 2 ); j >= 0; j-- )
    {
        r[j] -= gam[j + 1] * r[j + 1];
        //     cout << " j " << j << endl;
    }
}
#if ( PITTPACKACC )
#pragma acc routine seq
#endif
void TriDiag::enforceZeroMean( double *tmpMG, double *rh, double *diag, int index )
{

    int N = nChunk * nzChunk;
       
        rh[0]=0.0;
        rh[N-1]=0.0; 
        thomasLowMemNoBCV2( tmpMG, rh, diag, index );

}




#if ( PITTPACKACC )
#pragma acc routine seq
#endif
void TriDiag::shermanMorrisonThomas( double *tmpMG, double *rh, double *rh1, double diag, const double alpha, const double beta, int index )
{
    int N = nChunk * nzChunk;

    // print the input
    //
    /*
        cout<<" before solve "<<endl;
        for ( int i = 0; i < N; i++ )
        {
            cout<<rh[i]<<endl;
        }
    */

    double b[3];
    double bb[3];

// first enforce dirchlet type bc
// the goal is [1 ...... -1] at first and
// [-1 ...........1] at the last row
#if ( PITTPACKACC )
#pragma acc loop seq
#endif
    {
// this is due to the fact that -1 0 ..... 1 =0 in periodic
        b[0] = diag;
        b[1] = diag;
        b[2] = diag;
    }
    // cout<<diag<<endl;

    double fact, gamma;
    gamma = -b[0];

    // this is for diagonal entry modification
    // only 3 elements are needed
    bb[0] = b[0] - gamma;

    bb[1] = b[1];

    bb[2] = b[2] - alpha * beta / gamma;
    // enforcing boundary conditions here
    // note that supdiga[0]=0.0 and rh=0.0

/*    
        rh[0]=0.0;

        supDiag[0]=0.0;
    */

    // this will remove the singularity for the corner that we set the calue as zero
int counter=0;
#if(1)
    if ( fabs( diag + 2. ) < 1.e-10 )
    {
      //    cout<<" first solve's diag "<<diag<<endl;
      // enforceZeroMean( tmpMG, rh, b, index );      
      //  for ( int i = 0; i < N; i++ )
      //  {
      //          cout<<rh[i]<<endl;
      //  }
//cout<<RED<<" singlularity  "<<counter<<RESET<<endl;
// only when there is singularity remove it
      shermanMorrisonThomasV1( tmpMG, rh, rh1, diag, 1.0, 1.0, index );

      counter++;
  return;
    }
#endif
    thomasLowMemNoBC( tmpMG, rh, bb, index );
//  rhs is the new rhs
#if ( PITTPACKACC )
#pragma acc loop seq
#endif
    {
        rh1[0]     = gamma;
        rh1[N - 1] = alpha;
    }
    // no need for this, if already set to zero before calling this routine
    // inside solveThmBacth since this segment will not scale

#if ( PITTPACKACC )
#pragma acc loop seq
#endif
    for ( int i = 1; i < N - 1; i++ )
    {
        rh1[i] = 0.0;
        //    tmpMG[i]=0.0;
    }

    thomasLowMemNoBC( tmpMG, rh1, bb, index );
    
/*
       cout<<" second solve "<<endl;
        for ( int i = 0; i < N; i++ )
        {
            cout<<rh1[i]<<endl;
        }
  */  

    // move this for gang parallelism
    // in the future
    double part1 = ( rh[0] + beta * rh[N - 1] / gamma );
    double part2 = ( 1. + rh1[0] + beta * rh1[N - 1] / gamma );
/*
if(fabs(part2<1.e-6))
{
 // cout<<" index "<<index <<" "<<diag<<" "<<part2 <<endl;
}
*/
#if ( PITTPACKACC )
#pragma acc loop seq
#endif
    {
        // fact = ( rh[0] + beta * rh[N - 1] / gamma ) / ( 1. + rh1[0] + beta * rh1[N - 1] / gamma );
        fact = part1 / part2;
    }
#if ( 1 )
#if ( PITTPACKACC )
#pragma acc loop seq
#endif
    for ( int i = 0; i < N; i++ )
    {
        rh[i] -= fact * rh1[i];
        // rh[i] = rh[i]-  rh1[i]/part2;
        // rh[i] += rh1[i] ;
    }

#endif

    //    cout<<" ends index "<<index<<" eig "<<diag<<" "<<rh[0]<<" "<<rh[N-1]<<endl;
    //    cout<<" ends index "<<index<<" eig "<<diag<<" "<<endl;
/*    
        cout<<" final solve "<<endl;

        for ( int i = 0; i < N; i++ )
        {
            cout<<rh[i]<<endl;
        }
  */  
 #if(0)   
    if(fabs(rh[0]-rh[N-1])>1.e-6)
    {
     printf(" periodicity screwed up\n ");
     exit(0);
   }
#endif
  
}
#if ( PITTPACKACC )
#pragma acc routine seq
#endif
void TriDiag::shermanMorrisonThomasV1( double *tmpMG, double *rh, double *rh1, double diag, const double alpha, const double beta, int index )
{
    int N = nChunk * nzChunk;

    // print the input
    //
    /*
        cout<<" before solve "<<endl;
        for ( int i = 0; i < N; i++ )
        {
            cout<<rh[i]<<endl;
        }
    */

    double b[3];
    double bb[3];

// first enforce dirchlet type bc
// the goal is [1 ...... -1] at first and
// [-1 ...........1] at the last row
#if ( PITTPACKACC )
#pragma acc loop seq
#endif
    {
// this is due to the fact that -1 0 ..... 1 =0 in periodic
        b[0] = 1.0;
        b[1] = diag;
        b[2] = diag;
    }
    // cout<<diag<<endl;

    double fact, gamma;
    gamma = -b[0];

    // this is for diagonal entry modification
    // only 3 elements are needed
    bb[0] = b[0] - gamma;

    bb[1] = b[1];

    bb[2] = b[2] - alpha * beta / gamma;
    // enforcing boundary conditions here
    // note that supdiga[0]=0.0 and rh=0.0

   /* 
        rh[0]=0.0;
 
        supDiag[0]=0.0;
    */

    // this will remove the singularity for the corner that we set the calue as zero
#if(0)

    if ( fabs( diag + 2. ) < 1.e-10 )
    {
          cout<<"rhs before first solve's "<<diag<<endl;
        for ( int i = 0; i < N; i++ )
        {
                cout<<rh[i]<<endl;
        }
  //      return;
    }
#endif
    thomasLowMemNoBCV1( tmpMG, rh, bb, index );
//  rhs is the new rhs
#if ( PITTPACKACC )
#pragma acc loop seq
#endif
    {
        rh1[0]     = gamma;
        rh1[N - 1] = alpha;
    }
    // no need for this, if already set to zero before calling this routine
    // inside solveThmBacth since this segment will not scale

#if ( PITTPACKACC )
#pragma acc loop seq
#endif
    for ( int i = 1; i < N - 1; i++ )
    {
        rh1[i] = 0.0;
        //    tmpMG[i]=0.0;
    }

    thomasLowMemNoBCV1( tmpMG, rh1, bb, index );
    
/*
       cout<<" second solve "<<endl;
        for ( int i = 0; i < N; i++ )
        {
            cout<<rh1[i]<<endl;
        }
  */  

    // move this for gang parallelism
    // in the future
    double part1 = ( rh[0] + beta * rh[N - 1] / gamma );
    double part2 = ( 1. + rh1[0] + beta * rh1[N - 1] / gamma );
/*
if(fabs(part2<1.e-6))
{
 // cout<<" index "<<index <<" "<<diag<<" "<<part2 <<endl;
}
*/
#if ( PITTPACKACC )
#pragma acc loop seq
#endif
    {
        // fact = ( rh[0] + beta * rh[N - 1] / gamma ) / ( 1. + rh1[0] + beta * rh1[N - 1] / gamma );
        fact = part1 / part2;
    }
#if ( 1 )
#if ( PITTPACKACC )
#pragma acc loop seq
#endif
    for ( int i = 0; i < N; i++ )
    {
        rh[i] -= fact * rh1[i];
        // rh[i] = rh[i]-  rh1[i]/part2;
        // rh[i] += rh1[i] ;
    }

#endif

    //    cout<<" ends index "<<index<<" eig "<<diag<<" "<<rh[0]<<" "<<rh[N-1]<<endl;
    //    cout<<" ends index "<<index<<" eig "<<diag<<" "<<endl;
   #if(0) 
        cout<<" final solve "<<endl;

        for ( int i = 0; i < N; i++ )
        {
            cout<<rh[i]<<endl;
        } 
    if(fabs(rh[0]-rh[N-1])>1.e-6)
    {
     printf(" periodicity screwed up\n ");
     exit(0);
   }
#endif
  
}
