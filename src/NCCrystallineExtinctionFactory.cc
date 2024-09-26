
#include "NCCrystallineExtinctionFactory.hh"
#include "NCCrystallineExtinction.hh"

namespace NCPluginNamespace {

  class PluginScatter final : public NC::ProcImpl::ScatterIsotropicMat {
  public:

    //The factory wraps our custom CrystallineExtinction helper class in an NCrystal API
    //Scatter class.

    const char * name() const noexcept override { return NCPLUGIN_NAME_CSTR "Model"; }
    PluginScatter( CrystallineExtinction && pm ) : m_pm(std::move(pm)) {}

    NC::CrossSect crossSectionIsotropic(NC::CachePtr&, NC::NeutronEnergy ekin) const override
    {
      return NC::CrossSect{ m_pm.calcCrossSection( ekin.dbl() ) };
    }

    NC::ScatterOutcomeIsotropic sampleScatterIsotropic(NC::CachePtr&, NC::RNG& rng, NC::NeutronEnergy ekin ) const override
    {
      auto outcome = m_pm.sampleScatteringEvent( rng, ekin.dbl() );
      return { NC::NeutronEnergy{outcome.ekin_final}, NC::CosineScatAngle{outcome.mu} };
    }

  private:
    CrystallineExtinction m_pm;
  };

}

const char * NCP::CrystallineExtinctionFactory::name() const noexcept
{
  //Factory name. Keep this standardised form please:
  return NCPLUGIN_NAME_CSTR "Factory";
}

//////////////////////////////////////////////////////////////////////////////////
//                                                                              //
// Here follows the factory logic, for how the physics model provided by the    //
// plugin should be combined with existing models in NCrystal.                  //
//                                                                              //
// In the silly example here, we want our custom physics model to replace the   //
// existing incoherent-elastic model of NCrystal with our own model.            //
//                                                                              //
//////////////////////////////////////////////////////////////////////////////////

NC::Priority NCP::CrystallineExtinctionFactory::query( const NC::FactImpl::ScatterRequest& cfg ) const
{
  //Must return value >0 if we should do something, and a value higher than
  //100 means that we take precedence over the standard NCrystal factory:
  //if (!cfg.get_coh_elas())
  //  return NC::Priority::Unable;//coherent-elastic disabled, never do anything.

  //Ok, we might be applicable. Load input data and check if is something we
  //want to handle:

  if ( ! CrystallineExtinction::isApplicable( cfg.info() ) )
    return NC::Priority::Unable;

  //Ok, all good. Tell the framework that we want to deal with this, with a
  //higher priority than the standard factory gives (which is 100):
  return NC::Priority{999};
}

NC::ProcImpl::ProcPtr NCP::CrystallineExtinctionFactory::produce( const NC::FactImpl::ScatterRequest& cfg ) const
{
  //Ok, we are selected as the provider! First create our own scatter model:

  auto sc_ourmodel = NC::makeSO<PluginScatter>( CrystallineExtinction::createFromInfo( cfg.info() ) );

  //Now we just need to combine this with all the other physics
  //(i.e. Bragg+inelastic).  So ask the framework to set this up, except for
  //incoherent-elastic physics of course since we are now dealing with that
  //ourselves:

  auto sc_std = globalCreateScatter( cfg.modified("coh_elas=0") );

  //Combine and return:
  return combineProcs( sc_std, sc_ourmodel );
}
