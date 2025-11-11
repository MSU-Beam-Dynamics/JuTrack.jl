"""
Test script to verify pyJuTrack installation
"""

def test_import():
    """Test that pyJuTrack can be imported"""
    try:
        import pyJuTrack as jt
        print("✓ pyJuTrack imported successfully")
        return True
    except ImportError as e:
        print(f"✗ Failed to import pyJuTrack: {e}")
        return False

def test_basic_elements():
    """Test creating basic lattice elements"""
    try:
        import pyJuTrack as jt
        
        # Create elements
        d = jt.DRIFT("D1", 1.0)
        q = jt.KQUAD("Q1", length=0.5, k1=1.2)
        
        print("✓ Created DRIFT and KQUAD elements")
        return True
    except Exception as e:
        print(f"✗ Failed to create elements: {e}")
        return False

def test_lattice():
    """Test creating a lattice"""
    try:
        import pyJuTrack as jt
        
        d1 = jt.DRIFT("D1", 1.0)
        q1 = jt.KQUAD("Q1", length=0.5, k1=1.2)
        d2 = jt.DRIFT("D2", 1.0)
        
        lattice = jt.Lattice([d1, q1, d2])
        length = lattice.total_length()
        
        assert abs(length - 2.5) < 1e-10, f"Expected length 2.5, got {length}"
        
        print(f"✓ Created lattice with length {length:.2f} m")
        return True
    except Exception as e:
        print(f"✗ Failed to create lattice: {e}")
        return False

def test_beam():
    """Test creating a beam"""
    try:
        import numpy as np
        import pyJuTrack as jt
        
        # Create random particle distribution
        coords = np.random.randn(100, 6) * 1e-4
        beam = jt.Beam(coords, energy=1e9)
        
        print(f"✓ Created beam with {len(coords)} particles")
        return True
    except Exception as e:
        print(f"✗ Failed to create beam: {e}")
        return False

def test_tracking():
    """Test particle tracking"""
    try:
        import numpy as np
        import pyJuTrack as jt
        
        # Create simple lattice
        d1 = jt.DRIFT("D1", 1.0)
        q1 = jt.KQUAD("Q1", length=0.5, k1=1.2)
        d2 = jt.DRIFT("D2", 1.0)
        lattice = jt.Lattice([d1, q1, d2])
        
        # Create beam
        coords = np.random.randn(100, 6) * 1e-4
        beam = jt.Beam(coords, energy=1e9)
        
        # Track
        jt.linepass(lattice, beam)
        
        print("✓ Particle tracking successful")
        return True
    except Exception as e:
        print(f"✗ Tracking failed: {e}")
        return False

def test_tpsa():
    """Test TPSA functionality"""
    try:
        import pyJuTrack as jt
        
        # Set TPSA dimension
        jt.set_tps_dim(2)
        
        # Create DTPSAD variables
        x1 = jt.DTPSAD(1.0, 1)
        x2 = jt.DTPSAD(2.0, 2)
        
        # Simple computation
        result = x1 * x1 + x2 * x2
        
        print("✓ TPSA/DTPSAD functionality working")
        return True
    except Exception as e:
        print(f"✗ TPSA test failed: {e}")
        return False

def test_gradient():
    """Test automatic differentiation"""
    try:
        import pyJuTrack as jt
        
        jt.set_tps_dim(2)
        
        def objective(x1, x2):
            return x1 * x1 + x2 * x2 # python operator ** does not work with DTPSAD
        
        grad = jt.Gradient(objective, [1.0, 2.0])
        
        # Expected gradient: [2*x1, 2*x2] = [2.0, 4.0]
        assert abs(grad[0] - 2.0) < 1e-10, f"Expected grad[0]=2.0, got {grad[0]}"
        assert abs(grad[1] - 4.0) < 1e-10, f"Expected grad[1]=4.0, got {grad[1]}"
        
        print(f"✓ Automatic differentiation working (gradient: [{grad[0]:.1f}, {grad[1]:.1f}])")
        print("  Note: Python operators (+, -, *) work with DTPSAD!")
        print("  For power, use: jt.pow(x, n)")
        return True
    except Exception as e:
        print(f"✗ Gradient test failed: {e}")
        import traceback
        traceback.print_exc()
        return False

def run_all_tests():
    """Run all tests"""
    tests = [
        ("Import", test_import),
        ("Basic Elements", test_basic_elements),
        ("Lattice Creation", test_lattice),
        ("Beam Creation", test_beam),
        ("Particle Tracking", test_tracking),
        ("TPSA/DTPSAD", test_tpsa),
        ("Automatic Differentiation", test_gradient),
    ]
    
    results = []
    for name, test_func in tests:
        print(f"\nTest: {name}")
        print("-" * 70)
        try:
            success = test_func()
            results.append((name, success))
        except Exception as e:
            print(f"✗ Unexpected error: {e}")
            results.append((name, False))
    
    # Summary
    print("\n" + "=" * 70)
    print("Test Summary")
    print("=" * 70)
    
    passed = sum(1 for _, success in results if success)
    total = len(results)
    
    for name, success in results:
        status = "PASS" if success else "FAIL"
        symbol = "✓" if success else "✗"
        print(f"{symbol} {name:.<50} {status}")
    
    print("-" * 70)
    print(f"Total: {passed}/{total} tests passed")
    
    if passed == total:
        print("\n All tests passed! pyJuTrack is ready to use.")
        return 0
    else:
        print(f"\n {total - passed} test(s) failed. Check error messages above.")
        return 1

if __name__ == "__main__":
    import sys
    exit_code = run_all_tests()
    sys.exit(exit_code)
