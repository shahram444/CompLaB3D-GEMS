"""Shared XML fragments for the CompLB3D example suite.

Every example is assembled from these, so a change to the house style or a
correction to a tag name happens in one place instead of fifteen.

Note on XML comments: a comment body may not contain two consecutive hyphens.
That rules out the usual "-->" arrow and any "--" dash inside prose. Every
comment here is written accordingly; the validator checks it.
"""

NX, NY, NZ = 24, 24, 6


def header(number, title, lines):
    """The banner every case opens with."""
    body = "\n".join("     " + l if l else "" for l in lines)
    return """<?xml version="1.0" ?>
<!-- =========================================================================
     CompLB3D EXAMPLE %02d  |  %s
     =========================================================================

%s
     ========================================================================= -->
<parameters>
""" % (number, title, body)


def path(input_path="input"):
    return """
    <path>
        <src_path>src</src_path>
        <input_path>%s</input_path>
        <output_path>output</output_path>
    </path>
""" % input_path


def domain(geometry, materials="", peclet=0, delta_P="1e-6", tau="0.8",
           ade_max_iT=400, ns_max_iT1=20000, ade_update_interval=1,
           ns_update_interval=10):
    return """
    <LB_numerics>
        <domain>
            <nx>%d</nx>
            <ny>%d</ny>
            <nz>%d</nz>
            <dx>10</dx>
            <unit>um</unit>
            <characteristic_length>%d</characteristic_length>
            <filename>%s</filename>
            <material_numbers>
                <pore>2</pore>
                <solid>0</solid>
                <bounce_back>1</bounce_back>%s
            </material_numbers>
        </domain>

        <delta_P>%s</delta_P>
        <Peclet>%s</Peclet>
        <tau>%s</tau>
        <track_performance>false</track_performance>

        <iteration>
            <ns_max_iT1>%d</ns_max_iT1>
            <ns_max_iT2>2000</ns_max_iT2>
            <ns_converge_iT1>1e-6</ns_converge_iT1>
            <ns_converge_iT2>1e-4</ns_converge_iT2>
            <ns_update_interval>%d</ns_update_interval>
            <ade_update_interval>%d</ade_update_interval>
            <ade_max_iT>%d</ade_max_iT>
            <ade_converge_iT>1e-8</ade_converge_iT>
        </iteration>
    </LB_numerics>
""" % (NX, NY, NZ, NY, geometry, materials, delta_P, peclet, tau,
       ns_max_iT1, ns_update_interval, ade_update_interval, ade_max_iT)


def substrate(index, name, c0, D="5e-10", left=("Dirichlet", "0."),
              right=("Dirichlet", "0."), immobile=False, comment=None):
    c = ("                                    <!-- %s -->" % comment) if comment else ""
    imm = "\n            <immobile>true</immobile>" if immobile else ""
    return """
        <substrate%d>%s
            <name_of_substrates>%s</name_of_substrates>
            <initial_concentration>%s</initial_concentration>
            <substrate_diffusion_coefficients>
                <in_pore>%s</in_pore>
                <in_biofilm>%s</in_biofilm>
            </substrate_diffusion_coefficients>%s
            <left_boundary_type>%s</left_boundary_type>
            <left_boundary_condition>%s</left_boundary_condition>
            <right_boundary_type>%s</right_boundary_type>
            <right_boundary_condition>%s</right_boundary_condition>
        </substrate%d>
""" % (index, c, name, c0, D, D, imm, left[0], left[1], right[0], right[1], index)


def microbe(index, name, solver, reaction, b0, decay="0.0",
            visc=None, Dbio=None, Kc=None, Vmax=None, extra="", comment=None):
    """One microbe block.

    visc: <viscosity_ratio_in_biofilm>. It must be present for exactly the
    microbes that have a <microbeN> entry under material_numbers, and absent
    for the rest. The parser counts them and compares against that entry count,
    so an extra or a missing one terminates the run.

    Dbio: <biomass_diffusion_coefficients>. Mandatory when solver_type is FD.
    """
    c = ("                                <!-- %s -->" % comment) if comment else ""
    s = """
        <microbe%d>%s
            <name_of_microbes>%s</name_of_microbes>
            <solver_type>%s</solver_type>
            <reaction_type>%s</reaction_type>
            <initial_densities>%s</initial_densities>
            <decay_coefficient>%s</decay_coefficient>
""" % (index, c, name, solver, reaction, b0, decay)
    if visc is not None:
        s += "            <viscosity_ratio_in_biofilm>%s</viscosity_ratio_in_biofilm>\n" % visc
    if Dbio is not None:
        s += """            <biomass_diffusion_coefficients>
                <in_pore>%s</in_pore>
                <in_biofilm>%s</in_biofilm>
            </biomass_diffusion_coefficients>
""" % (Dbio, Dbio)
    if Kc is not None:
        s += "            <half_saturation_constants>%s</half_saturation_constants>\n" % Kc
    if Vmax is not None:
        s += "            <maximum_uptake_flux>%s</maximum_uptake_flux>\n" % Vmax
    s += extra
    s += """            <left_boundary_type>Neumann</left_boundary_type>
            <left_boundary_condition>0.</left_boundary_condition>
            <right_boundary_type>Neumann</right_boundary_type>
            <right_boundary_condition>0.</right_boundary_condition>
        </microbe%d>
""" % index
    return s


def io(vtk=100, chk=0):
    return """
    <IO>
        <read_NS_file>false</read_NS_file>
        <read_ADE_file>false</read_ADE_file>
        <ns_filename>nsLattice</ns_filename>
        <mask_filename>maskLattice</mask_filename>
        <subs_filename>subsLattice</subs_filename>
        <bio_filename>bioLattice</bio_filename>
        <save_VTK_interval>%d</save_VTK_interval>
        <save_CHK_interval>%d</save_CHK_interval>
        <debug_updRxn>0</debug_updRxn>
    </IO>

</parameters>
""" % (vtk, chk)


OFF_EQ = """
    <equilibrium>
        <enabled>false</enabled>
    </equilibrium>
"""

OFF_PRECIP = """
    <precipitation>
        <enabled>false</enabled>
    </precipitation>
"""

OFF_DISSOL = """
    <dissolution>
        <enabled>false</enabled>
    </dissolution>
"""
