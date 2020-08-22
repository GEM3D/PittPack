#include "definitions.h"
#include "pencilDcmp.hpp"

//#if ( PITTPACKACC )
void PoissonGPU::performTransformXdir() /*!< Is called on Host and Runs on GPU*/
{
#if ( PITTPACKACC )
    cufftHandle plan = NULL;
    double *    ptr  = P.P;

#pragma acc host_data use_device( ptr )
    {
        cufftPlan1d( &plan, nx, CUFFT_Z2Z, nyChunk * nzChunk );
        if ( CUFFT_STREAMS )
        {
            void *stream = acc_get_cuda_stream( acc_async_sync );
            cufftSetStream( plan, (cudaStream_t)stream );
        }
        if ( CUFFT_SUCCESS != cufftExecZ2Z( plan, (cufftDoubleComplex *)( ptr ), (cufftDoubleComplex *)( ptr ), CUFFT_FORWARD ) )
        {
            cout << " Exit Code : " << PittPackGetErrorEnum( CUFFT_FAIL_X ) << endl;
            exit( CUFFT_FAIL_X );
        }
        cufftDestroy( plan );
    }
#endif
}

void PoissonGPU::performInverseTransformXdir() /*!< Called on Host and Ran on GPU*/
{
#if ( PITTPACKACC )
    cufftHandle plan = NULL;
    double *    ptr  = P.P;

#pragma acc host_data use_device( ptr )
    {
        cufftPlan1d( &plan, nx, CUFFT_Z2Z, nyChunk * nzChunk );
        if ( CUFFT_STREAMS )
        {
            void *stream = acc_get_cuda_stream( acc_async_sync );
            cufftSetStream( plan, (cudaStream_t)stream );
        }
        if ( CUFFT_SUCCESS != cufftExecZ2Z( plan, (cufftDoubleComplex *)( ptr ), (cufftDoubleComplex *)( ptr ), CUFFT_INVERSE ) )
        {
            cout << " Exit Code : " << PittPackGetErrorEnum( CUFFT_FAIL_INV_X ) << endl;
            exit( CUFFT_FAIL_INV_X );
        }
        cufftDestroy( plan );
    }
#endif
}
void PoissonGPU::performTransformYdir() /*!< Called on Host and Ran on GPU*/
{
#if ( PITTPACKACC )
    cufftHandle plan = NULL;
    double *    ptr  = P.P;

#pragma acc host_data use_device( ptr )
    {
        cufftPlan1d( &plan, ny, CUFFT_Z2Z, nxChunk * nzChunk );
        if ( CUFFT_STREAMS )
        {
            void *stream = acc_get_cuda_stream( acc_async_sync );
            cufftSetStream( plan, (cudaStream_t)stream );
        }
        if ( CUFFT_SUCCESS != cufftExecZ2Z( plan, (cufftDoubleComplex *)( ptr ), (cufftDoubleComplex *)( ptr ), CUFFT_FORWARD ) )
        {
            cout << " Exit Code : " << PittPackGetErrorEnum( CUFFT_FAIL_Y ) << endl;
            exit( CUFFT_FAIL_Y );
        }
        cufftDestroy( plan );
    }
#endif
}

void PoissonGPU::performInverseTransformYdir() /*!< Called on Host and Ran on GPU*/
{
#if ( PITTPACKACC )
    cufftHandle plan = NULL;
    double *    ptr  = P.P;

#pragma acc host_data use_device( ptr )
    {
        cufftPlan1d( &plan, ny, CUFFT_Z2Z, nxChunk * nzChunk );
        if ( CUFFT_STREAMS )
        {
            void *stream = acc_get_cuda_stream( acc_async_sync );
            cufftSetStream( plan, (cudaStream_t)stream );
        }
        if ( CUFFT_SUCCESS != cufftExecZ2Z( plan, (cufftDoubleComplex *)( ptr ), (cufftDoubleComplex *)( ptr ), CUFFT_INVERSE ) )
        {
            cout << " Exit Code : " << PittPackGetErrorEnum( CUFFT_FAIL_INV_Y ) << endl;
            exit( CUFFT_FAIL_INV_Y );
        }
        cufftDestroy( plan );
    }
#endif
}

void PoissonGPU::triDiagCusparse( double *dl, double *ds, double *du, double *rhs )
{
    // cout<< " size of nz" <<nz<<endl;

#if ( PITTPACKACC )
    cusparseHandle_t handle = NULL;

#pragma acc host_data use_device( rhs, du, dl, ds )
    {
        cusparseCreate( &handle );
        if ( PIVOT == 0 )
        {
            cusparseDgtsv_nopivot( handle, nz, 1, dl, ds, du, rhs, nz );
        }
        else
        {
            cusparseDgtsv( handle, nz, 1, dl, ds, du, rhs, nz );
        }
    }

    cusparseDestroy( handle );
#endif
}


//#endif
