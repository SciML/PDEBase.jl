using PrecompileTools

@setup_workload begin
    # Minimal imports needed for precompilation
    @compile_workload begin
        # Create simple symbolic variables for precompilation
        @variables x y t
        @variables u(..)

        # Get unwrapped symbolic types
        x_sym = Symbolics.unwrap(x)
        y_sym = Symbolics.unwrap(y)
        t_sym = Symbolics.unwrap(t)
        u_term = Symbolics.unwrap(u(x, y, t))

        # Precompile safe_unwrap for common types
        safe_unwrap(u(x, y, t))
        safe_unwrap(x)
        safe_unwrap(1.0)

        # Create a differential term for testing
        Dx = Differential(x)
        term = Dx(u(x, y, t))
        term_sym = Symbolics.unwrap(term)

        # Precompile has_derivatives
        has_derivatives(term_sym)
        has_derivatives(u_term)
        has_derivatives(x_sym)

        # Precompile count_differentials
        count_differentials(term_sym, x_sym)
        count_differentials(u_term, x_sym)

        # Precompile differential_order
        differential_order(term_sym, x_sym)
        differential_order(u_term, x_sym)

        # Precompile get_depvars
        depvar_ops = [operation(u_term)]
        get_depvars(term_sym, depvar_ops)
        get_depvars(u_term, depvar_ops)

        # Precompile split_terms
        eq = Symbolics.unwrap(Dx(u(x, y, t)) ~ u(x, y, t))
        split_terms(eq)

        # Precompile split_additive_terms
        split_additive_terms(eq)

        # Precompile recursive_unwrap
        recursive_unwrap(term_sym)
    end
end
