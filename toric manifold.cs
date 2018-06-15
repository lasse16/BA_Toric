using System.Collections;
using System.Collections.Generic;
using UnityEngine;

public class toricmanifold : MonoBehaviour {

	private fixAngle _alpha;
	private fixAngle _theta;
	private fixAngle _phi;
	private GameObject _target1;
	private GameObject _target2;
	
	private Vector3 A = _target1.transform.position;
	private	Vector3 B = _target2.transform.position;
	private	Vector3 AB = A - B;
	
	public toricmanifold (float alpha, float theta, float phi, GameObject target1, GameObject target2) {
		_alpha = alpha;
		_phi = phi;
		_theta = theta;
		_target1 = target1;
		_target2 = target2;
	}
		public toricmanifold (float alpha, GameObject target1, GameObject target2) {
		_alpha = alpha;
		_phi = 0;
		_theta = 0;
		_target1 = target1;
		_target2 = target2;
	}
	
	public Vector3 ToWorldPosition(){
		Vector3 C = new Vector3(0,0,0);
		float last = Mathf.Sin(_alpha.toRad() + _theta.toRad()/2);
		
		Vector3 n = -AB;
		Vector2 n2 = n.projectZ();
		n2.rotate90();
		Vector3	z = Vector3(n2.x(),n2.y(),0); z.normalize();
		Vector3 t = Vector3.Cross(z,n);
		
		Quaternion qT = Quaternion.AngleAxis(t, _theta) ;
		Quaternion qP = Quaternion.AngleAxis(AB, _phi);
		C = A + (qP * qT * AB) * last;
		return C;
	} 
	
}
