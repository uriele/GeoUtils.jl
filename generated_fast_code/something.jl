         sinθ² = sinθ * sinθ
                  half_sin2θ = sinθ * cosθ
                  cos2θ = cosθ² - sinθ²
                end
                # compute the gradient
                begin
                  begin
                    dpx2dθ = -sinθ
                    dpy2dθ = bcosθ
                    ddx2dθ_0 = -bsinθ
                    ddy2dθ_0 = cosθ
                  end
                  begin
                    half_dR = e² * half_sin2θ
                    half_dR² = half_dR * half_dR
                    half_d²R = e² * cos2θ
                    N² = N * N
                    dNdθ = -half_dR
                    d²Ndθ² = -half_d²R + 3 * N² * half_dR²
                    ddx2dθ = (ddx2dθ_0 + N² * dNdθ * dx2) * N
                    ddy2dθ = (ddy2dθ_0 + N² * dNdθ * dy2) * N
                  end
                  begin
                    dFxdθ = dpx2dθ + s * ddx2dθ
                    dFydθ = dpy2dθ + s * ddy2dθ
                    dtdθ = dx1 * dFxdθ + dy1 * dFydθ
                    dPxdθ = dx1 * dtdθ
                    dPydθ = dy1 * dtdθ
                    dfxdθ = dFxdθ - dPxdθ
                    dfydθ = dFydθ - dPydθ
                  end
