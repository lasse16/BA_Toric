using System;
using System.Collections;
using System.Collections.Generic;
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

    private GameObject _target1;
    private GameObject _target2;
    private Camera _main;

    //target positions
    private Vector3 A;
    private Vector3 B;
    private Vector3 vecAB;


    //desired screen position 
    private Vector2 pA;
    private Vector2 pB;

    /**
     *Intialising all fields
     * @param alpha the angle between the camera and the two targets
     * @param theta horizontal rotation angle |gets automatically clamped between 1 and getMaxTheta()
     * @param phi   vertical rotation angle 
     * @param target1   first target
     * @param target2   second target
     */
    public Toricmanifold(float alpha, float theta, float phi, GameObject target1, GameObject target2)
    {
        _target1 = target1;
        _target2 = target2;
        A = _target1.transform.position;
        B = _target2.transform.position;
        vecAB = A - B;

        _alpha = new FixAngle(alpha);
        _phi = new FixAngle(phi);
        theta = Mathf.Clamp(theta, 1, getMaxTheta());
        _theta = new FixAngle(theta);

        //Better way to get the camera?
        _main = Camera.main;
        
    }

    /**
     * overloaded constructor with theta and phi = 0
     * 
     * 
     */

    public Toricmanifold(float alpha, GameObject target1, GameObject target2) : this(alpha, 0, 0, target1, target2)
    {
    }

    /**
    *Calculates the WorldPosition from the field values of the object
    * @return the vector3 position in world coordinates
    * 
    */
    public Vector3 ToWorldPosition()
    {
        Vector3 C;
        float last = (vecAB.magnitude * Mathf.Sin(_alpha.toRad() + _theta.toRad() / 2))/ Mathf.Sin(_alpha.toRad());


        Vector3 n = -vecAB;
       
        //n.projectZ
        Vector2 n2 = new Vector2(n.x, n.z);

        //n2.rotate90();
        float tmp = n2[0];
        n2[0] = -n2[1];
        n2[1] = tmp;

        Vector3 z = new Vector3(n2.x, 0, n2.y);
       
        
        Vector3 t = Vector3.Cross(z, n);

        Quaternion qT = Quaternion.AngleAxis(_theta.angle(), t);
        Quaternion qP = Quaternion.AngleAxis(_phi.angle(), n);
        C = A + (qP * qT * vecAB) * last / vecAB.magnitude;
        return C;
    }

    /**
     * Computes the orientation for a given position, optionally enforces a given tilt angle
     * @param camPos the camera position
     * @param TiltAngle = 0 the optional tilt angle around the cameras forward vector
     * 
     * @return Quaternion the final camera orientation
     * TODO check qTrans and qLook(mostly right)
     */
    public Quaternion ComputeOrientation(Vector3 camPos, float TiltAngle = 0)
    {
        //main camera direction
        Quaternion qLook = computeLookAt(camPos);

        //positions the targets at the desired screen positions
        Quaternion qTrans = computeLookAtTransition();

        //Tilt works
        Quaternion qPhi = computeTiltAngle(TiltAngle);

        return qPhi * qLook * (Quaternion.Inverse(qTrans));
    }

    /**
     * computes the rotation needed to place the targets at the desired screen positions
     * @return Quaternion with the  rotation
     * 
     * TODO check if pixel or procent for position specification //relative to camera (0.0) to (1,1)
     **/
    private Quaternion computeLookAtTransition()
    {
        
		//Lino Version
	
          Vector3  forward, up;
		{

            Vector3 pA3 = Vector3.Normalize(ToricComputing.GetVectorInCameraSpace(pA));
            Vector3 pB3 = Vector3.Normalize(ToricComputing.GetVectorInCameraSpace(pB));
			
			up = Vector3.Cross(pB3,pA3).normalized;		
			forward = (pA3 + pB3).normalized;
			
		}
		return Quaternion.LookRotation(forward,up);
         

        /**
         * 
         * 
        Vector2 pO = Vector2.zero;
        Vector2 pM = new Vector2((pA.x + pB.x) / 2, (pA.y + pB.y) / 2);
        Vector3 p3O = _main.ViewportToWorldPoint(Vector3.forward);
        Vector3.Normalize(p3O);

        //TODO correct way to calculate Vector3?
        Vector3 p3M = _main.ViewportToWorldPoint(new Vector3(pM.x / Sx, pM.y / Sy, 1));
        Vector3.Normalize(p3M);

        Vector3 m = Vector3.Cross(p3M, p3O);
        float angle = Vector3.Angle(p3M, p3O);

        return Quaternion.AngleAxis(angle, m);
        */
    }

    /**
    *enforces the tilt angle around the cameras forward vector
    * @param tilt the tilt angle
    * @return Quaternion the rotation aorund the forward vector
    */
    private Quaternion computeTiltAngle(float tilt)
    {
        FixAngle phi = new FixAngle(tilt);
        return Quaternion.AngleAxis(phi.angle(), _main.transform.forward);
    }

    /**
     * focuses on the middle between the two targets
     * @param camPos camera position
     * @return Quaternion the rotation to look at the middle between the two targets
     */
    private Quaternion computeLookAt(Vector3 camPos)
    {
        Vector3 dA = -1* (camPos - A);
        Vector3 dB = -1* (camPos - B);

       

        dA = Vector3.Normalize(dA);
        dB = Vector3.Normalize(dB);

        Vector3 l = 0.5f * (dA + dB);
        return Quaternion.LookRotation(l);

    }

    /**
     * sets the desired screen positions of each target
     * @param screenPos1 
     * @param screenPos2
     */
    public void SetDesiredPosition(Vector2 screenPos1, Vector2 screenPos2)
    {
        pA = screenPos1;
        pB = screenPos2;
    }

    //Returns the maximum value of theta for a given alpha
    public float getMaxTheta()
    {
        return 2 * (Mathf.PI - _alpha.toRad()) * Mathf.Rad2Deg;
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

    public GameObject getTarget1()
    {
        return _target1;
    }

    public GameObject getTarget2()
    {
        return _target2;
    }


    override
    public String ToString()
    {
        return "Alpha: " + _alpha.angle() + "Theta: " + _theta.angle() + "Phi: " + _phi.angle() + "Visibility: " + ToricComputing.visibilityCheck(this);
    }



    //TESTMETHODS


    public Quaternion testBasicLookAt(Vector3 campos)
    {
        return computeLookAt(campos);
    }


 




}

 