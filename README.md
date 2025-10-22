GENERAL INPUT GUIDELINES FOR CONDUCTION STUDIO

📋 GENERAL RULES (All Geometries)
Temperature Inputs:

✅ Inner temperature (T₁) should be greater than outer temperature (T₂) for heat to flow outward
✅ Use realistic temperature ranges (-200°C to 3000°C)
❌ Avoid equal temperatures (T₁ = T₂) → no heat transfer!
💡 Example: T₁ = 100°C, T₂ = 20°C ✓

Material Selection:

✅ Choose appropriate materials for your application
💡 Metals (Copper, Aluminum) for high conductivity
💡 Insulators (Foam, Fiberglass) for low conductivity
❌ Don't mix materials unrealistically (e.g., all insulation with high heat flux expectations)

Number of Layers:

✅ Use 1-5 layers (slider range)
💡 Start with 2 layers for composite problems
💡 Use 1 layer for simple single-material cases


🔲 PLANE WALL SPECIFIC RULES
Thickness (L):

✅ Must be positive (L > 0)
✅ Typical range: 0.001 m to 1.0 m
❌ Don't use zero or negative thickness
💡 Recommended: 0.01 m to 0.1 m for most problems
💡 Insulation layers: 0.02-0.05 m
💡 Metal layers: 0.005-0.05 m

Area (A):

✅ Default is 1.0 m² (built into code)
💡 Results scale linearly with area

Example Valid Setup:
Layer 1: Steel, L = 0.05 m
Layer 2: Insulation, L = 0.03 m
T₁ = 100°C, T₂ = 20°C
```

---

### **⭕ CYLINDER SPECIFIC RULES**

**Radii (r_inner, r_outer):**
- ✅ Must satisfy: **r_outer > r_inner > 0**
- ✅ For Layer 2, r_inner (Layer 2) = r_outer (Layer 1)
- ❌ Don't overlap or leave gaps between layers
- 💡 Typical pipe: r_inner = 0.02 m, r_outer = 0.025 m
- 💡 With insulation: add 0.03-0.05 m to outer radius

**Layer Continuity:**
- ✅ Each layer must connect: r₂(layer i) = r₁(layer i+1)
- ❌ No gaps: r_outer(L1)=0.03, r_inner(L2)=0.04 ❌
- ✅ Continuous: r_outer(L1)=0.03, r_inner(L2)=0.03 ✓

**Length (L):**
- ✅ Default is 1.0 m (built into code)
- 💡 Results scale linearly with length

**Example Valid Setup:**
```
Layer 1: Steel pipe
  r_inner = 0.02 m
  r_outer = 0.03 m
  
Layer 2: Insulation
  r_inner = 0.03 m  ← matches Layer 1 outer!
  r_outer = 0.06 m
  
T₁ = 100°C, T₂ = 20°C
```

**❌ Common Mistakes:**
```
Bad: Layer 1 (0.02 → 0.03), Layer 2 (0.04 → 0.05)  ← GAP!
Bad: Layer 1 (0.02 → 0.03), Layer 2 (0.025 → 0.04) ← OVERLAP!
```

---

### **🔵 SPHERE SPECIFIC RULES**

**Radii:**
- ✅ **r_outer must be > r_inner**
- ✅ For **solid sphere**: set r_inner = 0
- ✅ For **hollow sphere**: r_inner > 0 (e.g., 0.01 m)
- 💡 Typical range: r_outer = 0.01 m to 0.5 m

**Heat Generation (q̇):**
- ✅ Must be **positive** (q̇ > 0)
- ✅ Use scientific notation for large values: 1e6 = 1,000,000
- 💡 Typical range: 1×10⁴ to 1×10⁷ W/m³
- 💡 Nuclear fuel: ~1×10⁸ W/m³
- 💡 Electrical heating: ~1×10⁵ W/m³
- ❌ Don't use zero (no heat generation = no problem!)

**Convection Coefficient (h):**
- ✅ Must be **positive** (h > 0)
- 💡 Natural air convection: 5-25 W/m²·K
- 💡 Forced air convection: 25-250 W/m²·K
- 💡 Boiling water: 2500-100,000 W/m²·K
- 💡 Typical default: 100 W/m²·K

**Ambient Temperature (T∞):**
- ✅ Should be **less than** expected surface temperature
- 💡 Room temperature: 20-25°C
- 💡 For cooling problems, ensure T∞ < T_surface

**Example Valid Setup (Solid):**
```
Material: Copper
r_inner = 0 m          ← solid sphere
r_outer = 0.05 m
q̇ = 1e6 W/m³
h = 100 W/m²·K
T∞ = 25°C
```

**Example Valid Setup (Hollow):**
```
Material: Steel
r_inner = 0.01 m       ← hollow sphere
r_outer = 0.05 m
q̇ = 5e5 W/m³
h = 50 W/m²·K
T∞ = 20°C

🔗 CONTACT RESISTANCE RULES
Selection:

✅ Choose from the dropdown menu
✅ "Perfect Contact" (R″ = 0) for ideal cases
💡 Use realistic combinations (Steel/Steel, Aluminum/Aluminum)
💡 Higher pressure → lower contact resistance
❌ Don't select contact resistance for the last layer (no interface after it)

When to Use:

✅ Between dissimilar materials
✅ When surfaces are not perfectly bonded
✅ For realistic engineering problems
⚠️ Can significantly affect results!


⚠️ COMMON ERRORS & FIXES
Error MessageLikely CauseFix"could not convert string to float"Non-numeric inputEnter numbers only (use . for decimals)Negative heat fluxT₂ > T₁Swap temperaturesVery large/small resultsUnit mismatchCheck: meters (not mm), °C (not K)"math domain error"r_inner ≥ r_outerEnsure r_outer > r_inner

✅ QUICK VALIDATION CHECKLIST
Before clicking "Calculate", verify:
Plane Wall:

 All thicknesses > 0
 T₁ > T₂
 Reasonable material choices

Cylinder:

 r_outer > r_inner for each layer
 No gaps: r_out(layer i) = r_in(layer i+1)
 T₁ > T₂
 Innermost radius > 0

Sphere:

 r_outer > r_inner (or r_inner = 0 for solid)
 q̇ > 0
 h > 0
 Reasonable T∞ (< expected T_surface)


💡 BEST PRACTICES

Start Simple: Use 1-2 layers initially
Realistic Values: Stay within typical engineering ranges
Check Units: All inputs are in SI units (m, W, °C, W/m³)
Use Defaults: Start with default values and modify gradually
Compare Results: Do results make physical sense?

High k → low ΔT across layer
Low k (insulation) → high ΔT across layer
Heat flows from hot to cold


Experiment Safely:

Save working configurations
Use "Reset" button to return to defaults
Try "Show Steps" to understand calculations




📊 TYPICAL VALUE RANGES SUMMARY
ParameterMinimumTypicalMaximumTemperature-100°C20-500°C2000°CThickness (wall)0.001 m0.01-0.1 m1 mRadius0.001 m0.01-0.5 m5 mHeat generation1×10³ W/m³1×10⁵-10⁶ W/m³1×10⁸ W/m³Convection coeff.1 W/m²·K10-100 W/m²·K10,000 W/m²·K
