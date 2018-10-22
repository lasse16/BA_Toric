using System;
using UnityEngine;
/**
 * Representation of a toric manifold surface
 * Used for the main library funcctions like calculating the orientation, or the worldposition.
 * 
 * @author Lasse Haffke 
 * @version 20.06.18
 * 
 * 
 */
public class Toricmanifold
{
    //characteristics of toric manifold
    private FixAngle _alpha;
    private FixAngle _theta;
    private FixAngle _phi;

    public readonly GameObject _target1;
    public readonly GameObject _target2;
    private Camera _mainCamera;

    //target positions
    private Vector3 A;
    private Vector3 B;
    private Vector3 vecAB;


    //desired screen position 
    private Vector2 screenPositionA;
    private Vector2 screenPositionB;


    private Boolean screenPosSet;

    /// <summary>
    /// Sets up a toric manifold object with alpha, theta and phi.
    /// </summary>
    /// <remarks>
    /// alpha [ 1 - 180]°
    /// theta [ 1 - 2 * (180-alpha)]°
    ///  phi [ -180 - 180]°
    /// </remarks>
    /// <param name="alpha"> alpha</param>
    /// <param name="theta"> theta</param>
    /// <param name="phi"> phi</param>
    /// <param name="target1"> the first, left-most target</param>
    /// <param name="target2">the second target</param>
    public Toricmanifold(float alpha, float theta, float phi, GameObject target1, GameObject target2)
    {
        _target1 = target1;
        _target2 = target2;
        A = _target1.transform.position;
        B = _target2.transform.position;
        vecAB = B - A;

        _alpha = new FixAngle(alpha);
        _phi = new FixAngle(phi);
        theta = Mathf.Clamp(theta, 1, getMaxTheta());
        _theta = new FixAngle(theta);

        _mainCamera = Camera.main;
    }

    public Toricmanifold(float alpha, GameObject target1, GameObject target2) : this(alpha, 0, 0, target1, target2)
    {
    }

    /// <summary>
    /// Converts the toric representation to world coordinates
    /// </summary>
    /// <returns>Vector3 in world coordinates</returns>
    public Vector3 ToWorldPosition()
    {
        Vector3 C;
        float last = (vecAB.magnitude * Mathf.Sin(_alpha.toRad() + _theta.toRad() / 2)) / Mathf.Sin(_alpha.toRad());

        Vector3 n = -vecAB;
        n = n.normalized;

        Vector2 n2 = new Vector2(n.x, n.z);

        float tmp = n2[0];
        n2[0] = -n2[1];
        n2[1] = tmp;

        Vector3 z = new Vector3(n2.x, 0, n2.y);
        z = z.normalized;

        Vector3 t = Vector3.Cross(z, n);

        //horizontal rotation
        Quaternion qT = Quaternion.AngleAxis(_theta.angle() / 2, t.normalized);

        //vertical rotation
        Quaternion qP = Quaternion.AngleAxis(_phi.angle(), n);

        Vector3 res = qP * qT * vecAB;
        C = res * last / vecAB.magnitude + A;
        return C;
    }


    

    /// <summary>
    /// Computes the proper orientation to put the targets at the screen positions
    /// </summary>
    /// <param name="camPos"></param>
    /// <param name="TiltAngle">user-wished tilt angle in Degrees</param>
    /// <returns></returns>
    public Quaternion ComputeOrientation(Vector3 camPos, float TiltAngle = 0)
    {
        if (!screenPosSet) throw new Exception("Please set the screen positions");

        Quaternion qLook = computeLookAt(camPos);

        Quaternion qTrans = computeLookAtTransition();

        Quaternion qPhi = computeTiltAngle(TiltAngle);

        return qPhi * qLook * (Quaternion.Inverse(qTrans));
    }




    /// <summary>
    /// Sets the desired screen positions
    /// </summary>
    /// <param name="screenPos1">on-screen position target 1</param>
    /// <param name="screenPos2">on-screen position target 1</param>
    public void SetDesiredPosition(Vector2 screenPos1, Vector2 screenPos2)
    {
        screenPosSet = true;
        screenPositionA = screenPos1;
        screenPositionB = screenPos2;
    }

    //Returns the maximum value of theta for a given alpha
    public float getMaxTheta()
    {
        return 2 * (180 - _alpha.angle());
    }

    public float getAlpha()
    {
        return _alpha.angle();
    }

    public float getTheta()
    {
        return _theta.angle();
    }

    public float getPhi()
    {
        return _phi.angle();
    }

    override
    public String ToString()
    {
        return "Alpha: " + _alpha.angle() + "Theta: " + _theta.angle() + "Phi: " + _phi.angle() + "Visibility: " + ToricComputing.visibilityCheck(this);
    }



    //private methods

    private Quaternion computeLookAtTransition()
    {
        Vector3 forward, up;
        {

            Vector3 pA3 = Vector3.Normalize(ToricComputing.GetVectorInCameraSpace(screenPositionA));
            Vector3 pB3 = Vector3.Normalize(ToricComputing.GetVectorInCameraSpace(screenPositionB));

            up = Vector3.Cross(pB3, pA3).normalized;
            forward = (pA3 + pB3).normalized;

        }
        return Quaternion.LookRotation(forward, up);
    }


    private Quaternion computeTiltAngle(float tilt)
    {
        FixAngle phi = new FixAngle(tilt);
        return Quaternion.AngleAxis(phi.angle(), _mainCamera.transform.forward);
    }

    private Quaternion computeLookAt(Vector3 camPos)
    {
        Vector3 dA = (A - camPos).normalized;
        Vector3 dB = (B - camPos).normalized;

        Vector3 l = 0.5f * (dA + dB);
        return Quaternion.LookRotation(l);

    }
}

