
    struct InitGradPerTaskData
    {
      unsigned int                              d;
      unsigned int                              vel_dpc;
      unsigned int                              pres_dpc;
      FullMatrix<double>                        local_grad;
      std::vector<types::global_dof_index>      vel_local_dof_indices;
      std::vector<types::global_dof_index>      pres_local_dof_indices;
      InitGradPerTaskData (const unsigned int   dd,
                           const unsigned int   vdpc,
                           const unsigned int   pdpc )
        :
        d(dd),
        vel_dpc (vdpc),
        pres_dpc (pdpc),
        local_grad (vdpc, pdpc),
        vel_local_dof_indices (vdpc),
        pres_local_dof_indices (pdpc)
      {}
    };


    struct AdvectionScratchData
    {
      unsigned int                              nqp;
      unsigned int                              dpc;
      std::vector< Point<dim> >                 u_star_local;
      std::vector< Tensor<1,dim> >              grad_u_star;
      std::vector<double>                       u_star_tmp;
      FEValues<dim>                             fe_val;
      AdvectionScratchData (const FE_Q<dim> &   fe,
                            const QGauss<dim> & quad,
                            const UpdateFlags   flags )
        :
        nqp (quad.size()),
        dpc (fe.dofs_per_cell),
        u_star_local (nqp),
        grad_u_star (nqp),
        u_star_tmp (nqp),
        fe_val (fe, quad, flags)
      {}
      AdvectionScratchData (const AdvectionScratchData &data)
        :
        nqp (data.nqp),
        dpc (data.dpc),
        u_star_local (nqp),
        grad_u_star (nqp),
        u_star_tmp (nqp),
        fe_val (data.fe_val.get_fe(),
                data.fe_val.get_quadrature(),
                data.fe_val.get_update_flags())
      {}
    };

    struct InitGradScratchData
    {
      unsigned int  nqp;
      FEValues<dim> fe_val_vel;
      FEValues<dim> fe_val_pres;
      InitGradScratchData (const FE_Q<dim> &    fe_v,
                           const FE_Q<dim> &    fe_p,
                           const QGauss<dim> &  quad,
                           const UpdateFlags    flags_v,
                           const UpdateFlags    flags_p )
        :
        nqp (quad.size()),
        fe_val_vel (fe_v, quad, flags_v),
        fe_val_pres (fe_p, quad, flags_p)
      {}
      InitGradScratchData (const InitGradScratchData &data)
        :
        nqp (data.nqp),
        fe_val_vel (    data.fe_val_vel.get_fe(),
                        data.fe_val_vel.get_quadrature(),
                        data.fe_val_vel.get_update_flags()),
        fe_val_pres (   data.fe_val_pres.get_fe(),
                        data.fe_val_pres.get_quadrature(),
                        data.fe_val_pres.get_update_flags())
      {}
    };


    struct AdvectionPerTaskData
    {
      FullMatrix<double>                        local_advection;
      std::vector<types::global_dof_index>      local_dof_indices;
      AdvectionPerTaskData (const unsigned int  dpc     )
        :
        local_advection (dpc, dpc),
        local_dof_indices (dpc)
      {}
    };